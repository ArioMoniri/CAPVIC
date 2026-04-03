"""FastMCP server with 26 variant classification tools.

Multi-source precision oncology and variant classification MCP server.
Integrates CIViC, ClinVar, OncoKB, VICC MetaKB, gnomAD, UniProt, PubMed, and MyVariant.info.
"""

from __future__ import annotations

import asyncio
import logging
import sys
from typing import Any

from mcp.server.fastmcp import FastMCP

from variant_mcp.classification import ACMGAMPHelper, AMPTierClassifier, OncogenicityScorer
from variant_mcp.clients.base_client import ClientError
from variant_mcp.clients.civic_client import CIViCClient
from variant_mcp.clients.clinvar_client import ClinVarClient
from variant_mcp.clients.gnomad_client import GnomADClient
from variant_mcp.clients.metakb_client import MetaKBClient
from variant_mcp.clients.myvariant_client import MyVariantClient
from variant_mcp.clients.oncokb_client import ONCOKB_NO_TOKEN_MSG, OncoKBClient
from variant_mcp.clients.pubmed_client import PubMedClient
from variant_mcp.clients.uniprot_client import UniProtClient
from variant_mcp.constants import DISCLAIMER
from variant_mcp.formatters.reports import ReportFormatter
from variant_mcp.formatters.tables import TableFormatter
from variant_mcp.models.evidence import EvidenceBundle
from variant_mcp.utils.variant_normalizer import VariantNormalizer

logger = logging.getLogger("variant_mcp")
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter("%(asctime)s [%(name)s] %(levelname)s: %(message)s"))
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# --- Clients (module-level singletons) ---
civic_client = CIViCClient()
clinvar_client = ClinVarClient()
oncokb_client = OncoKBClient()
metakb_client = MetaKBClient()
gnomad_client = GnomADClient()
pubmed_client = PubMedClient()
uniprot_client = UniProtClient()
myvariant_client = MyVariantClient()

# --- Utilities ---
variant_normalizer = VariantNormalizer()

# --- Classification engines ---
amp_classifier = AMPTierClassifier()
oncogenicity_scorer = OncogenicityScorer()
acmg_helper = ACMGAMPHelper()

# --- Formatters ---
report_fmt = ReportFormatter()
table_fmt = TableFormatter()

# --- MCP Server ---
mcp = FastMCP(
    "variant_mcp",
    instructions=(
        "Multi-source precision oncology and variant classification MCP server.\n\n"
        "Integrates CIViC, ClinVar, OncoKB, and VICC MetaKB to provide comprehensive "
        "clinical evidence for genetic variant interpretation. Implements AMP/ASCO/CAP "
        "somatic tiering, ClinGen/CGC/VICC oncogenicity scoring, and ACMG/AMP germline "
        "pathogenicity reference.\n\n"
        "KEY TOOLS:\n"
        "- variant_search_evidence: Multi-source evidence search (start here)\n"
        "- variant_classify: Apply classification frameworks (somatic or germline)\n"
        "- variant_pathogenicity_summary: Pathogenic vs Benign evidence comparison\n"
        "- score_oncogenicity: ClinGen/CGC/VICC point-based scoring\n"
        "- classify_amp_tier: AMP/ASCO/CAP somatic tiering\n"
        "- civic_search_evidence / clinvar_search / oncokb_annotate: Source-specific\n"
        "- lookup_gnomad_frequency: Population allele frequencies from gnomAD\n"
        "- lookup_protein_domains: Protein domain mapping from UniProt\n"
        "- predict_variant_effect: In-silico predictions (SIFT, PolyPhen-2, REVEL, etc.)\n"
        "- search_literature / get_publication: PubMed literature search\n"
        "- normalize_variant: HGVS notation parsing and conversion\n"
        "- lookup_gene / lookup_disease / lookup_therapy: Discovery/autocomplete\n"
        "- get_classification_frameworks_reference: Understand the guidelines\n\n"
        "IMPORTANT: All outputs are for research/education only, not clinical decisions.\n"
        "OncoKB features require a free academic API token (ONCOKB_API_TOKEN env var)."
    ),
)


# ============================================================================
# Helper: gather evidence from all sources
# ============================================================================
async def _gather_evidence(
    gene: str,
    variant: str | None = None,
    disease: str | None = None,
    sources: list[str] | None = None,
    limit: int = 25,
) -> EvidenceBundle:
    """Query all available data sources and return an EvidenceBundle."""
    bundle = EvidenceBundle(gene=gene, variant=variant or "", disease=disease)
    query_all = sources is None
    tasks: dict[str, Any] = {}

    if query_all or "civic" in (sources or []):
        tasks["civic"] = civic_client.search_evidence_parsed(
            gene=gene, variant=variant, disease=disease, first=limit
        )
    if query_all or "clinvar" in (sources or []):
        tasks["clinvar"] = clinvar_client.search_variants(
            gene=gene, variant_name=variant, disease=disease, limit=limit
        )
    if (query_all or "oncokb" in (sources or [])) and oncokb_client.is_available and variant:
        tasks["oncokb"] = oncokb_client.annotate_mutation(gene, variant)
    if query_all or "metakb" in (sources or []):
        tasks["metakb"] = metakb_client.search(gene=gene, variant=variant, disease=disease)

    # Phase 10: UniProt domain + in-silico predictions (when variant is provided)
    if variant and (query_all or "uniprot" in (sources or [])):
        tasks["uniprot"] = uniprot_client.check_variant_in_domain(gene, variant)
    if variant and (query_all or "myvariant" in (sources or [])):
        tasks["myvariant"] = myvariant_client.search_variant(gene, variant)

    results = await asyncio.gather(
        *[tasks[k] for k in tasks],
        return_exceptions=True,
    )

    for key, result in zip(tasks.keys(), results, strict=False):
        if isinstance(result, Exception):
            bundle.errors[key] = f"{type(result).__name__}: {result}"
            logger.error("Error from %s: %s", key, result)
        elif key == "civic":
            bundle.civic_evidence = result  # type: ignore[assignment]
        elif key == "clinvar":
            bundle.clinvar_variants = result  # type: ignore[assignment]
        elif key == "oncokb":
            if isinstance(result, str):
                bundle.errors["oncokb"] = result
            else:
                bundle.oncokb_annotation = result  # type: ignore[assignment]
        elif key == "metakb":
            bundle.metakb_interpretations = result  # type: ignore[assignment]
        elif key == "uniprot":
            if hasattr(result, "domains"):
                bundle.protein_domains = result.domains  # type: ignore[union-attr]
        elif key == "myvariant":
            if result is not None:
                bundle.in_silico_predictions = result  # type: ignore[assignment]

    # Also fetch CIViC assertions for the gene
    if query_all or "civic" in (sources or []):
        try:
            bundle.civic_assertions = await civic_client.search_assertions(
                gene=gene, disease=disease, first=10
            )
        except Exception as e:
            bundle.errors["civic_assertions"] = str(e)

    return bundle


# ============================================================================
# CATEGORY 1: Multi-Source Search Tools
# ============================================================================


@mcp.tool(
    annotations={  # type: ignore[arg-type]
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": True,
    }
)
async def variant_search_evidence(
    gene: str,
    variant: str | None = None,
    disease: str | None = None,
    therapy: str | None = None,
    evidence_type: str | None = None,
    sources: list[str] | None = None,
    limit: int = 25,
) -> str:
    """Perform a comprehensive multi-database search for clinical evidence on genetic variants, therapies, and diseases.

    This is the PRIMARY search tool — queries CIViC, ClinVar, OncoKB, MetaKB, UniProt, and MyVariant.info simultaneously, then returns a unified evidence report with summary statistics, top evidence listings, clinical interpretations, and citation sources.

    Use this for questions like:
    - "Find evidence for colorectal cancer therapies involving KRAS mutations"
    - "What is the clinical significance of BRAF V600E in melanoma?"
    - "Are there any FDA-approved therapies targeting EGFR T790M?"

    The returned report includes structured markdown tables and evidence breakdowns suitable for data visualization and comparison.

    Args:
        gene: Gene symbol (e.g., BRAF, KRAS, TP53, EGFR). Required.
        variant: Variant name (e.g., V600E, G12C, T790M). Optional — omit to get all variants for a gene.
        disease: Disease or cancer type (e.g., "colorectal cancer", "melanoma", "NSCLC"). Optional.
        therapy: Therapy or drug name filter (e.g., "vemurafenib", "sotorasib"). Optional.
        evidence_type: Filter by evidence type: PREDICTIVE (therapy response), DIAGNOSTIC, PROGNOSTIC, ONCOGENIC, FUNCTIONAL. Optional.
        sources: List of databases to query: civic, clinvar, oncokb, metakb, uniprot, myvariant. Default: all.
        limit: Max results per source (1-100, default 25).
    """
    bundle = await _gather_evidence(gene, variant, disease, sources, limit)
    return report_fmt.format_evidence_report(bundle)


@mcp.tool(
    annotations={  # type: ignore[arg-type]
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": True,
    }
)
async def variant_classify(
    gene: str,
    variant: str,
    variant_origin: str,
    disease: str | None = None,
) -> str:
    """Apply formal variant classification frameworks to produce a structured clinical assessment with evidence codes, confidence levels, and tier assignments.

    For somatic variants: applies both AMP/ASCO/CAP 4-tier classification (Li et al. 2017) AND ClinGen/CGC/VICC oncogenicity SOP point-based scoring (Horak et al. 2022). Returns tier assignment (I-IV), evidence level (A-D), oncogenicity score, and all applied evidence codes with point values.

    For germline variants: retrieves ACMG/AMP 5-tier interpretation from ClinVar with review status and submitter breakdown.

    The output includes structured markdown tables ideal for reports and visualizations.

    Args:
        gene: Gene symbol (e.g., BRAF, TP53, KRAS). Required.
        variant: Variant name (e.g., V600E, R175H, G12C). Required.
        variant_origin: Must be "somatic" or "germline". Required.
        disease: Disease or cancer type context (e.g., "melanoma", "lung adenocarcinoma"). Optional but recommended for somatic variants.
    """
    from variant_mcp.models.classification import VariantClassificationReport

    origin = variant_origin.lower().strip()
    if origin not in ("somatic", "germline"):
        return f"Error: variant_origin must be 'somatic' or 'germline', got '{variant_origin}'"

    bundle = await _gather_evidence(gene, variant, disease)

    report = VariantClassificationReport(
        gene=gene,
        variant=variant,
        disease=disease,
        variant_origin=origin,
        sources_queried=bundle.sources_queried,
    )

    if origin == "somatic":
        report.amp_tier = amp_classifier.classify(bundle)
        report.oncogenicity = oncogenicity_scorer.score_variant(
            gene, variant, evidence_bundle=bundle
        )
    else:
        # Germline: ACMG/AMP interpretation from ClinVar
        if bundle.clinvar_variants:
            report.acmg_interpretation = acmg_helper.explain_clinvar_classification(
                bundle.clinvar_variants[0]
            )

    result = report_fmt.format_classification_report(report)

    # Suggest fetching gnomAD data if not in bundle (improves SBVS1/SBS1/OM4/OP4)
    if not bundle.has_gnomad_data:
        result += (
            "\n\n---\n\n> **Tip**: For more accurate oncogenicity scoring "
            "(SBVS1/SBS1/OM4/OP4 evidence codes), first fetch population frequency "
            "data using `lookup_gnomad_frequency` with the genomic coordinate "
            "(e.g., `7-140753336-A-T` for BRAF V600E on GRCh38), then re-run "
            "classification. The scorer uses gnomAD allele frequencies to apply "
            "frequency-based evidence codes per Horak et al. 2022."
        )

    return result


@mcp.tool(
    annotations={  # type: ignore[arg-type]
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": True,
    }
)
async def variant_compare_sources(
    gene: str,
    variant: str,
) -> str:
    """Cross-reference a variant across ALL databases and produce a side-by-side concordance/discordance analysis.

    Returns a structured comparison table showing how CIViC, ClinVar, OncoKB, and MetaKB classify the same variant, highlighting agreements and conflicts. Ideal for generating comparison visualizations and identifying data gaps between sources.

    Args:
        gene: Gene symbol (e.g., BRAF, EGFR). Required.
        variant: Variant name (e.g., V600E, L858R). Required.
    """
    bundle = await _gather_evidence(gene, variant)
    return report_fmt.format_source_comparison(bundle)


# ============================================================================
# CATEGORY 2: CIViC-Specific Tools
# ============================================================================


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def civic_search_evidence(
    gene: str | None = None,
    variant: str | None = None,
    disease: str | None = None,
    therapy: str | None = None,
    evidence_type: str | None = None,
    significance: str | None = None,
    limit: int = 25,
) -> str:
    """Search the CIViC (Clinical Interpretation of Variants in Cancer) database for curated clinical evidence linking genetic variants to therapies, diagnoses, and prognoses.

    Returns individual evidence items with evidence level (A-E), evidence type, clinical significance, associated therapies, diseases, and PubMed citations. Supports flexible filtering by gene, variant, disease, therapy, and evidence type.

    Use this for questions like "What therapies target KRAS G12C?" or "Find prognostic evidence for TP53 mutations in breast cancer."

    Args:
        gene: Gene symbol (e.g., BRAF, KRAS). Optional.
        variant: Variant name (e.g., V600E, G12C). Optional.
        disease: Disease name (e.g., "colorectal cancer", "melanoma"). Optional.
        therapy: Therapy or drug name (e.g., "sotorasib", "pembrolizumab"). Optional.
        evidence_type: Filter by type: PREDICTIVE (therapy response), DIAGNOSTIC, PROGNOSTIC, PREDISPOSING, ONCOGENIC, FUNCTIONAL. Optional.
        significance: Clinical significance filter (e.g., "Sensitivity", "Resistance"). Optional.
        limit: Max results (1-100, default 25).
    """
    try:
        items = await civic_client.search_evidence_parsed(
            gene=gene,
            variant=variant,
            disease=disease,
            therapy=therapy,
            evidence_type=evidence_type,
            first=limit,
        )
    except ClientError as e:
        return f"CIViC search error: {e}\n\n{DISCLAIMER}"

    if not items:
        parts = []
        if gene:
            parts.append(f"gene='{gene}'")
        if variant:
            parts.append(f"variant='{variant}'")
        if disease:
            parts.append(f"disease='{disease}'")
        return (
            f"CIViC returned 0 results for {', '.join(parts)}. "
            "Try broadening your search — remove variant/disease filters, "
            f"or check the gene symbol with lookup_gene.\n\n{DISCLAIMER}"
        )

    lines = [f"# CIViC Evidence Search Results ({len(items)} items)\n"]
    for item in items:
        therapies_str = ", ".join(item.therapies) if item.therapies else "N/A"
        pmid = f"PMID:{item.source.pmid}" if item.source and item.source.pmid else ""
        lines.append(
            f"- **EID{item.id}** [{item.evidence_type}] Level {item.evidence_level} — "
            f"{item.significance or 'N/A'} | Disease: {item.disease or 'N/A'} | "
            f"Therapies: {therapies_str} {pmid}"
        )
    lines.append(f"\n{DISCLAIMER}")
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def civic_search_assertions(
    gene: str | None = None,
    disease: str | None = None,
    therapy: str | None = None,
    significance: str | None = None,
    limit: int = 25,
) -> str:
    """Search CIViC curated assertions — expert-reviewed, higher-confidence clinical interpretations.

    Assertions represent the highest confidence tier in CIViC: each is reviewed by domain experts and includes AMP/ASCO/CAP tier assignments. Use this to find FDA-guideline-level evidence for specific gene/disease/therapy combinations.

    Use for questions like:
    - "Are there AMP Tier I assertions for BRAF in melanoma?"
    - "What therapies have curated assertions for EGFR-mutant NSCLC?"
    - "Show me predictive assertions for KRAS in colorectal cancer"

    Args:
        gene: Gene symbol (e.g., BRAF, KRAS, EGFR). Optional.
        disease: Disease name (e.g., "melanoma", "lung adenocarcinoma"). Optional.
        therapy: Therapy name (e.g., "vemurafenib", "pembrolizumab"). Optional.
        significance: Significance filter (e.g., "sensitivityresponse", "resistance"). Optional.
        limit: Max results (1-100, default 25).
    """
    try:
        assertions = await civic_client.search_assertions(
            gene=gene,
            disease=disease,
            therapy=therapy,
            significance=significance,
            first=limit,
        )
    except ClientError as e:
        return f"CIViC assertions error: {e}\n\n{DISCLAIMER}"

    if not assertions:
        return f"No CIViC assertions found for the given filters.\n\n{DISCLAIMER}"

    lines = [f"# CIViC Assertions ({len(assertions)} results)\n"]
    for a in assertions:
        lines.append(
            f"- **AID{a.id}** {a.assertion_type}: {a.significance} | "
            f"AMP: {a.amp_level or 'N/A'} | Disease: {a.disease or 'N/A'} | "
            f"Therapies: {', '.join(a.therapies) if a.therapies else 'N/A'} | "
            f"Evidence items: {a.evidence_count}"
        )
    lines.append(f"\n{DISCLAIMER}")
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def civic_get_gene(name: str) -> str:
    """Get gene details and all associated variants from the CIViC knowledge base.

    Returns the gene's official name, Entrez ID, description, and a list of all CIViC-curated variants with their IDs. Use this to discover which variants have clinical evidence for a given gene.

    Args:
        name: Gene symbol (e.g., BRAF, KRAS, EGFR, TP53).
    """
    try:
        data = await civic_client.get_gene(name)
    except ClientError as e:
        return f"CIViC gene lookup error: {e}\n\n{DISCLAIMER}"

    if not data:
        return f"Gene '{name}' not found in CIViC. Check spelling with lookup_gene.\n\n{DISCLAIMER}"

    variants = data.get("variants", {})
    variant_nodes = variants.get("nodes", [])
    total = variants.get("totalCount", 0)

    lines = [
        f"# CIViC Gene: {data.get('name', name)}\n",
        f"**Official Name**: {data.get('officialName', 'N/A')}",
        f"**Entrez ID**: {data.get('entrezId', 'N/A')}",
        f"**Total Variants**: {total}",
    ]
    if data.get("description"):
        lines.append(f"\n**Description**: {data['description'][:500]}")
    if variant_nodes:
        lines.append("\n**Variants**:")
        for v in variant_nodes[:20]:
            lines.append(f"- {v.get('name', 'Unknown')} (ID: {v.get('id')})")
    lines.append(f"\n{DISCLAIMER}")
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def civic_get_variant(variant_id: int) -> str:
    """Get detailed variant information from CIViC by variant ID.

    Returns the variant name, gene, variant types, and molecular profile score. Use after civic_search_evidence or civic_get_gene to drill into a specific variant.

    Args:
        variant_id: CIViC variant ID number (from search results or gene lookups).
    """
    try:
        data = await civic_client.get_variant(variant_id)
    except ClientError as e:
        return f"CIViC variant lookup error: {e}\n\n{DISCLAIMER}"

    if not data:
        return f"Variant ID {variant_id} not found in CIViC.\n\n{DISCLAIMER}"

    gene = data.get("gene", {})
    mp = data.get("singleVariantMolecularProfile", {}) or {}
    vtypes = data.get("variantTypes", [])

    lines = [
        f"# CIViC Variant: {data.get('name', 'Unknown')}\n",
        f"**Gene**: {gene.get('name', 'N/A')}",
        f"**Variant Types**: {', '.join(t.get('name', '') for t in vtypes) if vtypes else 'N/A'}",
        f"**Molecular Profile Score**: {mp.get('molecularProfileScore', 'N/A')}",
        f"\n{DISCLAIMER}",
    ]
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def civic_get_evidence_item(evidence_id: int) -> str:
    """Get full details for a single CIViC evidence item including clinical interpretation, source citation, and evidence rating.

    Returns evidence type, level, significance, direction, disease, therapies, rating (1-5 stars), description, and source citation with PMID. Use after civic_search_evidence to get the full text and context of a specific evidence item.

    Args:
        evidence_id: CIViC evidence item ID number (from search results).
    """
    try:
        item = await civic_client.get_evidence_item(evidence_id)
    except ClientError as e:
        return f"CIViC evidence item error: {e}\n\n{DISCLAIMER}"

    if not item:
        return f"Evidence item {evidence_id} not found in CIViC.\n\n{DISCLAIMER}"

    lines = [
        f"# CIViC Evidence Item EID{item.id}\n",
        f"**Type**: {item.evidence_type} | **Level**: {item.evidence_level} | "
        f"**Significance**: {item.significance}",
        f"**Gene**: {item.gene or 'N/A'} | **Variant**: {item.variant or 'N/A'}",
        f"**Disease**: {item.disease or 'N/A'}",
        f"**Therapies**: {', '.join(item.therapies) if item.therapies else 'N/A'}",
        f"**Direction**: {item.evidence_direction or 'N/A'}",
        f"**Rating**: {item.evidence_rating or 'N/A'}/5",
    ]
    if item.description:
        lines.append(f"\n**Description**: {item.description[:1000]}")
    if item.source:
        lines.append(f"\n**Citation**: {item.source.citation or 'N/A'}")
        if item.source.pmid:
            lines.append(f"**PMID**: {item.source.pmid}")
    lines.append(f"\n{DISCLAIMER}")
    return "\n".join(lines)


# ============================================================================
# CATEGORY 3: ClinVar Tools
# ============================================================================


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def clinvar_search(
    gene: str | None = None,
    variant: str | None = None,
    clinical_significance: str | None = None,
    disease: str | None = None,
    limit: int = 25,
) -> str:
    """Search NCBI ClinVar for germline and somatic variant classifications, expert review status, star ratings, and submitter consensus data.

    ClinVar aggregates interpretations from clinical laboratories worldwide. Returns classifications (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign), review stars (0-4), associated conditions, and protein changes. Now includes germline, oncogenicity, and clinical impact classifications.

    Args:
        gene: Gene symbol (e.g., BRCA1, TP53). Optional but recommended.
        variant: Variant name (e.g., C61G). Optional.
        clinical_significance: Filter by classification: pathogenic, likely pathogenic, uncertain significance, likely benign, benign. Optional.
        disease: Disease or phenotype filter (e.g., "breast cancer"). Optional.
        limit: Max results (1-100, default 25).
    """
    if not any([gene, variant, clinical_significance, disease]):
        return "Error: At least one search parameter is required.\n\n" + DISCLAIMER

    try:
        variants = await clinvar_client.search_variants(
            gene=gene,
            variant_name=variant,
            clinical_significance=clinical_significance,
            disease=disease,
            limit=limit,
        )
    except ClientError as e:
        return f"ClinVar search error: {e}\n\n{DISCLAIMER}"

    if not variants:
        return (
            f"ClinVar returned 0 results. Try broadening your search or checking "
            f"gene/variant names.\n\n{DISCLAIMER}"
        )

    lines = [f"# ClinVar Search Results ({len(variants)} variants)\n"]
    for cv in variants:
        lines.append(f"### {cv.title or 'Unknown'}")
        lines.append(f"- **Variation ID**: {cv.variation_id}")
        lines.append(
            f"- **Classification**: {cv.clinical_significance} "
            f"({cv.review_stars} — {cv.review_status})"
        )
        if cv.genes:
            lines.append(f"- **Genes**: {', '.join(cv.genes)}")
        if cv.conditions:
            lines.append(f"- **Conditions**: {', '.join(cv.conditions[:5])}")
        if cv.protein_change:
            lines.append(f"- **Protein Change**: {cv.protein_change}")
        if cv.clinvar_url:
            lines.append(f"- **URL**: {cv.clinvar_url}")
        lines.append("")
    lines.append(DISCLAIMER)
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def clinvar_get_variant(
    variation_id: int | None = None,
    rsid: str | None = None,
    hgvs: str | None = None,
) -> str:
    """Get a full ClinVar variant record with all submitter classifications, HGVS expressions, conditions, and conflict detection.

    Returns the aggregate classification, review status with star rating, all submitter interpretations with breakdown, HGVS expressions, associated conditions, and flags conflicting interpretations. Provide exactly one identifier.

    Use for questions like:
    - "What is the ClinVar classification for variation 65533?"
    - "Look up rs113488022 in ClinVar"
    - "Get ClinVar details for NM_004333.6:c.1799T>A"

    Args:
        variation_id: ClinVar variation ID number (e.g., 65533). Optional.
        rsid: dbSNP rsID (e.g., rs113488022). Optional.
        hgvs: HGVS expression (e.g., NM_004333.6:c.1799T>A). Optional.
    """
    if not any([variation_id, rsid, hgvs]):
        return "Error: Provide variation_id, rsid, or hgvs.\n\n" + DISCLAIMER

    try:
        if variation_id:
            cv = await clinvar_client.get_variant_detail(variation_id)
            variants = [cv] if cv else []
        elif rsid:
            variants = await clinvar_client.get_variant_by_rsid(rsid)
        else:
            assert hgvs is not None
            variants = await clinvar_client.get_variant_by_hgvs(hgvs)
    except ClientError as e:
        return f"ClinVar lookup error: {e}\n\n{DISCLAIMER}"

    if not variants:
        return f"No ClinVar record found.\n\n{DISCLAIMER}"

    cv = variants[0]
    lines = [f"# ClinVar Variant: {cv.title or 'Unknown'}\n"]
    lines.append(f"**Variation ID**: {cv.variation_id}")
    lines.append(
        f"**Classification**: {cv.clinical_significance} ({cv.review_stars} — {cv.review_status})"
    )
    if cv.genes:
        lines.append(f"**Genes**: {', '.join(cv.genes)}")
    if cv.rsid:
        lines.append(f"**rsID**: {cv.rsid}")
    if cv.conditions:
        lines.append(f"**Conditions**: {', '.join(cv.conditions)}")
    if cv.hgvs_expressions:
        lines.append(f"**HGVS**: {', '.join(cv.hgvs_expressions[:5])}")
    if cv.last_evaluated:
        lines.append(f"**Last Evaluated**: {cv.last_evaluated}")

    if cv.submitter_classifications:
        lines.append(f"\n**Submitter Classifications** ({cv.submitter_count} total):")
        breakdown: dict[str, int] = {}
        for sub in cv.submitter_classifications:
            cls = sub.get("classification", "Unknown")
            breakdown[cls] = breakdown.get(cls, 0) + 1
        for cls, count in sorted(breakdown.items(), key=lambda x: -x[1]):
            lines.append(f"- {cls}: {count} submitter(s)")

    if cv.conflicting:
        lines.append("\n**WARNING: CONFLICTING interpretations among submitters**")

    if cv.clinvar_url:
        lines.append(f"\n**URL**: {cv.clinvar_url}")
    lines.append(f"\n{DISCLAIMER}")
    return "\n".join(lines)


# ============================================================================
# CATEGORY 4: OncoKB Tools
# ============================================================================


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def oncokb_annotate(
    gene: str,
    variant: str,
    tumor_type: str | None = None,
) -> str:
    """Annotate a somatic mutation with OncoKB — the FDA-recognized precision oncology knowledge base.

    Returns oncogenicity classification, mutation effect (gain/loss-of-function), FDA therapeutic levels (1-4, R1-R2), specific drug recommendations, and gene/variant summaries. Requires a free academic API token (ONCOKB_API_TOKEN env var).

    Use for questions like "Is BRAF V600E actionable in melanoma?" or "What drugs target EGFR L858R?"

    Args:
        gene: Gene symbol (e.g., BRAF, EGFR, ALK). Required.
        variant: Protein change (e.g., V600E, L858R, T790M). Required.
        tumor_type: OncoTree tumor type code (e.g., MEL for melanoma, NSCLC, CRC). Optional but recommended — enables tumor-specific therapeutic levels.
    """
    if not oncokb_client.is_available:
        return ONCOKB_NO_TOKEN_MSG + "\n\n" + DISCLAIMER

    try:
        result = await oncokb_client.annotate_mutation(gene, variant, tumor_type)
    except ClientError as e:
        return f"OncoKB annotation error: {e}\n\n{DISCLAIMER}"

    if isinstance(result, str):
        return result + "\n\n" + DISCLAIMER

    from variant_mcp.constants import ONCOKB_LEVELS

    lines = [f"# OncoKB Annotation: {gene} {variant}\n"]
    lines.append(f"**Oncogenic**: {result.oncogenic or 'N/A'}")
    lines.append(f"**Mutation Effect**: {result.known_effect or 'N/A'}")
    if result.mutation_effect:
        lines.append(f"**Effect Description**: {result.mutation_effect[:500]}")
    if result.highest_sensitive_level:
        desc = ONCOKB_LEVELS.get(result.highest_sensitive_level, "")
        lines.append(f"**Highest Sensitive Level**: {result.highest_sensitive_level} ({desc})")
    if result.highest_resistance_level:
        desc = ONCOKB_LEVELS.get(result.highest_resistance_level, "")
        lines.append(f"**Highest Resistance Level**: {result.highest_resistance_level} ({desc})")
    if result.treatments:
        lines.append("\n**Treatments**:")
        for tx in result.treatments:
            drugs = ", ".join(tx.get("drugs", []))
            lines.append(f"- {tx.get('level', '')}: {drugs}")
    if result.gene_summary:
        lines.append(f"\n**Gene Summary**: {result.gene_summary[:500]}")
    if result.variant_summary:
        lines.append(f"**Variant Summary**: {result.variant_summary[:500]}")
    if result.oncokb_url:
        lines.append(f"\n**URL**: {result.oncokb_url}")
    lines.append(f"\n{DISCLAIMER}")
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def oncokb_cancer_genes() -> str:
    """Get the complete OncoKB curated cancer gene list with oncogene/tumor suppressor classification and therapeutic levels.

    Returns a structured table of all cancer genes with oncogene status, tumor suppressor gene (TSG) status, and highest sensitive therapeutic level. Use to check if a gene is a known cancer driver, or to compare multiple genes. Output is a markdown table ideal for visualizations. Requires ONCOKB_API_TOKEN environment variable.
    """
    if not oncokb_client.is_available:
        return ONCOKB_NO_TOKEN_MSG + "\n\n" + DISCLAIMER

    try:
        result = await oncokb_client.get_cancer_gene_list()
    except ClientError as e:
        return f"OncoKB error: {e}\n\n{DISCLAIMER}"

    if isinstance(result, str):
        return result + "\n\n" + DISCLAIMER

    lines = [f"# OncoKB Cancer Gene List ({len(result)} genes)\n"]
    lines.append("| Gene | Oncogene | TSG | Highest Level |")
    lines.append("|------|----------|-----|---------------|")
    for g in result[:100]:
        hugo = g.get("hugoSymbol", "")
        oncogene = "Yes" if g.get("oncogene") else "No"
        tsg = "Yes" if g.get("tsg") else "No"
        level = g.get("highestSensitiveLevel", "N/A")
        lines.append(f"| {hugo} | {oncogene} | {tsg} | {level} |")
    if len(result) > 100:
        lines.append(f"\n*Showing first 100 of {len(result)} genes.*")
    lines.append(f"\n{DISCLAIMER}")
    return "\n".join(lines)


# ============================================================================
# CATEGORY 5: Classification Framework Tools
# ============================================================================


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def classify_amp_tier(
    gene: str,
    variant: str,
    disease: str | None = None,
) -> str:
    """Apply the AMP/ASCO/CAP 4-tier somatic variant classification framework (Li et al. 2017) to assign clinical significance tiers and evidence levels.

    Gathers evidence from CIViC, ClinVar, and OncoKB, then assigns Tier I (strong significance, FDA-approved), Tier II (potential significance), Tier III (unknown), or Tier IV (benign). Returns the tier, evidence level (A-D), confidence assessment, evidence trail, and sources used.

    Args:
        gene: Gene symbol (e.g., BRAF, KRAS). Required.
        variant: Variant name (e.g., V600E, G12C). Required.
        disease: Disease context (e.g., "melanoma"). Optional but recommended — affects tier assignment.
    """
    bundle = await _gather_evidence(gene, variant, disease)
    result = amp_classifier.classify(bundle)
    return report_fmt.format_amp_tier(result) + f"\n\n{DISCLAIMER}"


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def score_oncogenicity(
    gene: str,
    variant: str,
    evidence_codes: list[str] | None = None,
) -> str:
    """Apply the ClinGen/CGC/VICC Oncogenicity SOP (Horak et al. 2022) point-based scoring system with 18 evidence codes.

    Scores a somatic variant on a scale from Benign (≤-7) through VUS (0-5) to Oncogenic (≥10). Auto-detects evidence codes from CIViC, ClinVar, OncoKB, gnomAD frequencies, UniProt domains, and in-silico predictors — or accepts manually specified codes. Returns total points, classification, individual evidence codes with explanations, and confidence level.

    Args:
        gene: Gene symbol (e.g., BRAF, TP53). Required.
        variant: Variant name (e.g., V600E, R175H). Required.
        evidence_codes: Explicit evidence codes (e.g., ['OVS1', 'OS3', 'OM4']). Optional — if omitted, auto-detects from available data.
    """
    bundle = None
    if not evidence_codes:
        bundle = await _gather_evidence(gene, variant)

    result = oncogenicity_scorer.score_variant(
        gene,
        variant,
        evidence_codes=evidence_codes,
        evidence_bundle=bundle,
    )
    return report_fmt.format_oncogenicity(result) + f"\n\n{DISCLAIMER}"


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def explain_acmg_criteria(
    criteria_code: str | None = None,
    query: str | None = None,
) -> str:
    """Reference tool for ACMG/AMP germline pathogenicity classification criteria (Richards et al. 2015).

    Look up any of the 28 ACMG/AMP criteria codes with full descriptions and strength levels. Provide a specific code for details, search by keyword, or use 'all' for the complete reference table.

    Use for questions like:
    - "What does PVS1 mean in ACMG classification?"
    - "Explain the PM2 criterion"
    - "Show me all ACMG pathogenicity criteria"

    Args:
        criteria_code: ACMG criteria code (e.g., PVS1, PM2, BS1, PP3). Optional.
        query: Keyword search or 'all' for complete reference. Optional.
    """
    if criteria_code:
        return acmg_helper.get_criteria_reference(criteria_code) + f"\n\n{DISCLAIMER}"
    if query and query.lower().strip() == "all":
        return acmg_helper.get_all_criteria() + f"\n\n{DISCLAIMER}"
    if query:
        # Try to find matching criteria
        from variant_mcp.constants import ACMG_CRITERIA

        matches = []
        for code, desc in ACMG_CRITERIA.items():
            if query.lower() in desc.lower() or query.lower() in code.lower():
                matches.append(f"**{code}**: {desc}")
        if matches:
            return "\n".join(
                ["# ACMG/AMP Criteria Matching Query\n"] + matches + [f"\n{DISCLAIMER}"]
            )
        return f"No ACMG criteria matching '{query}'. Use criteria_code='PVS1' or query='all'.\n\n{DISCLAIMER}"

    return (
        "Provide criteria_code (e.g., 'PVS1') or query='all' for the complete reference.\n\n"
        + DISCLAIMER
    )


# ============================================================================
# CATEGORY 6: Discovery Tools
# ============================================================================


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def lookup_gene(query: str) -> str:
    """Gene name autocomplete and discovery across the CIViC knowledge base.

    Type-ahead search that matches partial gene names. Use this to find the correct gene symbol before running evidence searches, or to explore which genes have clinical evidence in CIViC.

    Args:
        query: Gene name or partial name (e.g., 'BRA' matches BRAF, 'EGF' matches EGFR).
    """
    try:
        results = await civic_client.search_typeahead("gene", query)
    except ClientError as e:
        return f"Gene lookup error: {e}"

    if not results:
        return f"No genes matching '{query}' found in CIViC."

    lines = [f"# Gene Lookup: '{query}' ({len(results)} matches)\n"]
    for g in results:
        lines.append(f"- **{g.get('name', '')}** (Entrez: {g.get('entrezId', 'N/A')})")
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def lookup_disease(query: str) -> str:
    """Disease name autocomplete and discovery across the CIViC knowledge base.

    Type-ahead search that matches partial disease names and returns Disease Ontology IDs (DOID). Use to find the correct disease name before running evidence or classification searches.

    Args:
        query: Disease name or partial name (e.g., 'melan' matches Melanoma, 'colorect' matches Colorectal Cancer).
    """
    try:
        results = await civic_client.search_typeahead("disease", query)
    except ClientError as e:
        return f"Disease lookup error: {e}"

    if not results:
        return f"No diseases matching '{query}' found in CIViC."

    lines = [f"# Disease Lookup: '{query}' ({len(results)} matches)\n"]
    for d in results:
        doid = d.get("doid", "N/A")
        lines.append(f"- **{d.get('name', '')}** (DOID: {doid})")
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def lookup_therapy(query: str) -> str:
    """Therapy and drug name autocomplete and discovery across the CIViC knowledge base.

    Type-ahead search that matches partial therapy/drug names and returns NCI Thesaurus IDs. Use to find the correct therapy name before searching for treatment evidence, or to explore available targeted therapies.

    Args:
        query: Drug or therapy name or partial (e.g., 'vemur' matches Vemurafenib, 'pembroli' matches Pembrolizumab).
    """
    try:
        results = await civic_client.search_typeahead("therapy", query)
    except ClientError as e:
        return f"Therapy lookup error: {e}"

    if not results:
        return f"No therapies matching '{query}' found in CIViC."

    lines = [f"# Therapy Lookup: '{query}' ({len(results)} matches)\n"]
    for t in results:
        ncit = t.get("ncitId", "N/A")
        lines.append(f"- **{t.get('name', '')}** (NCIt: {ncit})")
    return "\n".join(lines)


# ============================================================================
# CATEGORY 7: Utility Tools
# ============================================================================


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def get_classification_frameworks_reference() -> str:
    """Get a comprehensive reference document explaining all three variant classification frameworks used in precision oncology and genetics.

    Returns a detailed guide covering AMP/ASCO/CAP 4-tier somatic classification (Li et al. 2017), ClinGen/CGC/VICC Oncogenicity SOP point-based scoring (Horak et al. 2022), and ACMG/AMP 5-tier germline pathogenicity (Richards et al. 2015). Explains when to use each framework, how they relate, and decision criteria. Ideal for understanding the classification methodology before interpreting results.
    """
    return table_fmt.format_frameworks_reference()


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def variant_pathogenicity_summary(
    gene: str,
    variant: str,
    disease: str | None = None,
) -> str:
    """Generate a structured pathogenic-vs-benign evidence report showing all supporting and opposing evidence for a variant's clinical significance.

    Aggregates data from CIViC, ClinVar, OncoKB, gnomAD, UniProt domains, and in-silico predictors into a clear two-column breakdown: evidence supporting pathogenicity/oncogenicity vs evidence supporting benign classification. Includes oncogenicity score and clinical interpretation. Output is structured markdown ideal for reports, artifacts, and visualizations.

    Args:
        gene: Gene symbol (e.g., BRAF, KRAS). Required.
        variant: Variant name (e.g., V600E, G12C). Required.
        disease: Disease context (e.g., "melanoma"). Optional.
    """
    bundle = await _gather_evidence(gene, variant, disease)
    oncogenicity = oncogenicity_scorer.score_variant(gene, variant, evidence_bundle=bundle)
    return report_fmt.format_pathogenicity_summary(bundle, oncogenicity)


# ============================================================================
# Phase 10: Scientific Enhancement Tools
# ============================================================================


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def lookup_gnomad_frequency(
    variant_id: str,
    genome_version: str = "GRCh38",
) -> str:
    """Look up population allele frequencies from the Genome Aggregation Database (gnomAD) — the largest human exome/genome frequency resource (807k exomes, 76k genomes).

    Returns global allele frequency, per-population breakdown (African, European, East Asian, South Asian, etc.), homozygote counts, and clinical interpretation. Critical for ACMG BA1 (>5% = benign standalone), PM2 (absent from controls), and oncogenicity SOP codes SBVS1/SBS1/OM4/OP4. Output includes a structured frequency table ideal for visualization.

    Args:
        variant_id: gnomAD variant ID in chrom-pos-ref-alt format.
            Use GRCh38 coordinates with GRCh38 (e.g., "7-140753336-A-T" for BRAF V600E).
            Use GRCh37 coordinates with GRCh37 (e.g., "7-140453136-A-T" for BRAF V600E).
        genome_version: "GRCh38" (default, gnomAD v4) or "GRCh37" (gnomAD v2.1).
    """
    try:
        freq = await gnomad_client.get_variant_frequency(variant_id, genome_version)

        lines = [f"## 🧬 gnomAD Population Frequency — {freq.variant_id or 'N/A'}\n"]
        if freq.rsid:
            lines.append(f"**rsID**: {freq.rsid}")
        lines.append(f"**Genome version**: {freq.genome_version}")
        lines.append(f"**Source**: {freq.source}")
        lines.append(f"**Filter status**: {freq.filter_status or 'N/A'}\n")
        af_str = f"{freq.allele_frequency:.6g}" if freq.allele_frequency is not None else "N/A"
        lines.append(f"**Global allele frequency**: {af_str}")
        lines.append(f"**Allele count**: {freq.allele_count or 'N/A'}")
        lines.append(f"**Allele number**: {freq.allele_number or 'N/A'}")
        lines.append(f"**Homozygote count**: {freq.homozygote_count or 'N/A'}\n")

        if freq.population_frequencies:
            lines.append("### Per-Population Frequencies\n")
            lines.append("| Population | Allele Frequency |")
            lines.append("|-----------|-----------------|")
            for pop, pop_af in sorted(freq.population_frequencies.items()):
                lines.append(f"| {pop} | {pop_af:.6g} |")

        af: float | None = freq.allele_frequency
        lines.append("\n### Clinical Interpretation\n")
        if af is not None:
            if af > 0.05:
                lines.append(
                    "🟢 **AF > 5%** → Supports **BA1** (ACMG benign standalone) "
                    "and **SBVS1** (oncogenicity SOP, −8 points)"
                )
            elif af > 0.01:
                lines.append("🔵 **AF > 1%** → Supports **SBS1** (oncogenicity SOP, −4 points)")
            elif af > 0.0001:
                lines.append(
                    "🟡 **AF 0.01–1%** → Does not meet benign frequency thresholds; "
                    "VUS-supporting range"
                )
            else:
                lines.append("🟠 **AF < 0.01%** → Supports **OM4** (rare in population, +2 points)")
        else:
            lines.append(
                "🔴 **Absent from gnomAD** → Supports **OP4** (+1 point) and "
                "**PM2** (ACMG absent from controls)"
            )

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)
    except ClientError as exc:
        return f"gnomAD query error: {exc}\n\n{DISCLAIMER}"


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def normalize_variant(
    variant: str,
    gene: str | None = None,
) -> str:
    """Parse and normalize variant notation between multiple HGVS and protein change formats.

    Converts between shorthand (V600E), 1-letter protein (p.V600E), 3-letter HGVS protein (p.Val600Glu), and cDNA notation. Detects variant type (missense, nonsense, frameshift, splice, in-frame deletion, etc.) and extracts position, reference, and alternate amino acids.

    Use for questions like:
    - "Convert V600E to HGVS notation"
    - "What type of variant is p.Arg175His?"
    - "Normalize the notation for c.1799T>A"

    Args:
        variant: Variant notation to parse (e.g., V600E, p.Val600Glu, c.1799T>A, R175H).
        gene: Gene symbol for context (optional, helps with gene-specific formatting).
    """
    notation = variant_normalizer.normalize(variant)

    lines = ["## 🔀 Variant Notation\n"]
    lines.append(f"**Input**: `{notation.original}`")
    if gene:
        lines.append(f"**Gene**: {gene}")
    lines.append(f"**Variant type**: {notation.variant_type or 'Unknown'}")
    if notation.position:
        lines.append(f"**Position**: {notation.position}")
    if notation.ref_aa:
        lines.append(f"**Reference AA**: {notation.ref_aa}")
    if notation.alt_aa:
        lines.append(f"**Alternate AA**: {notation.alt_aa}")

    lines.append("\n### Notation Formats\n")
    lines.append("| Format | Notation |")
    lines.append("|--------|----------|")
    if notation.protein_1letter:
        lines.append(f"| 1-letter protein | `{notation.protein_1letter}` |")
    if notation.protein_3letter:
        lines.append(f"| 3-letter HGVS protein | `{notation.protein_3letter}` |")
    if notation.cdna:
        lines.append(f"| cDNA | `{notation.cdna}` |")

    lines.append(f"\n---\n\n{DISCLAIMER}")
    return "\n".join(lines)


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def lookup_protein_domains(
    gene: str,
    variant: str | None = None,
    position: int | None = None,
) -> str:
    """Look up protein functional domains from UniProt/InterPro and check whether a variant falls within a critical functional region.

    Returns all annotated protein domains (kinase, DNA-binding, active sites, transmembrane, etc.), maps the variant position to overlapping domains, and assesses OM1 evidence (critical domain without benign variation). Output includes a structured domain table ideal for protein domain visualizations.

    Args:
        gene: Gene symbol. Required.
        variant: Variant name — position will be extracted automatically.
        position: Explicit amino acid position to check.
    """
    try:
        if variant:
            result = await uniprot_client.check_variant_in_domain(gene, variant)
        elif position:
            domains = await uniprot_client.get_domain_at_position(gene, position)
            from variant_mcp.models.evidence import DomainCheckResult

            result = DomainCheckResult(
                gene=gene,
                variant=f"position {position}",
                position=position,
                in_domain=bool(domains),
                domains=domains,
                evidence_for_om1=bool(domains),
            )
        else:
            features = await uniprot_client.get_protein_features(gene)
            lines = [f"## 🏗️ Protein Domains — {gene}\n"]
            lines.append(f"**UniProt ID**: {features.uniprot_id or 'N/A'}")
            lines.append(f"**Protein length**: {features.protein_length or 'N/A'} aa\n")
            if features.domains:
                lines.append("| Domain | Type | Start | End | Source |")
                lines.append("|--------|------|-------|-----|--------|")
                for d in features.domains:
                    lines.append(
                        f"| {d.name} | {d.domain_type or '—'} | {d.start_pos} | "
                        f"{d.end_pos} | {d.source or '—'} |"
                    )
            else:
                lines.append("No domains found in UniProt for this gene.")
            lines.append(f"\n---\n\n{DISCLAIMER}")
            return "\n".join(lines)

        lines = [f"## 🏗️ Domain Check — {gene} {result.variant}\n"]
        lines.append(f"**Position**: {result.position or 'N/A'}")
        lines.append(f"**In functional domain**: {'✅ Yes' if result.in_domain else '❌ No'}")
        lines.append(f"**Supports OM1**: {'✅ Yes' if result.evidence_for_om1 else '❌ No'}\n")

        if result.domains:
            lines.append("### Overlapping Domains\n")
            lines.append("| Domain | Type | Range | Source |")
            lines.append("|--------|------|-------|--------|")
            for d in result.domains:
                lines.append(
                    f"| {d.name} | {d.domain_type or '—'} | "
                    f"{d.start_pos}–{d.end_pos} | {d.source or '—'} |"
                )
            lines.append(
                "\n**OM1**: Located in critical/functional domain without benign "
                "variation (+2 points in oncogenicity SOP)"
            )
        else:
            lines.append("No known functional domains overlap this position.")

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)
    except ClientError as exc:
        return f"UniProt query error: {exc}\n\n{DISCLAIMER}"


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def search_literature(
    gene: str,
    variant: str | None = None,
    disease: str | None = None,
    limit: int = 10,
) -> str:
    """Search PubMed for peer-reviewed publications about a gene, variant, or disease combination.

    Returns total publication count and recent papers with titles, authors, journals, years, PMIDs, and DOIs. Use to find supporting literature for variant interpretation, review the evidence base for a gene-disease association, or compile a bibliography.

    Use for questions like:
    - "Find recent publications about BRAF V600E in melanoma"
    - "How many papers exist on KRAS mutations in colorectal cancer?"
    - "Search for literature on TP53 and Li-Fraumeni syndrome"

    Args:
        gene: Gene symbol (e.g., BRAF, KRAS, TP53). Required.
        variant: Variant name to narrow search (e.g., V600E). Optional.
        disease: Disease context to narrow search (e.g., "melanoma"). Optional.
        limit: Max publications to return (1-50, default 10).
    """
    try:
        result = await pubmed_client.search_publications(gene, variant, disease, limit)

        lines = ["## 📚 PubMed Literature Search\n"]
        lines.append(f"**Query**: `{result.query}`")
        lines.append(f"**Total publications found**: {result.total_count}\n")

        if result.publications:
            for i, pub in enumerate(result.publications, 1):
                authors = ", ".join(pub.authors[:3])
                if len(pub.authors) > 3:
                    authors += " et al."
                lines.append(f"### {i}. {pub.title or 'Untitled'}\n")
                lines.append(f"**Authors**: {authors}")
                lines.append(f"**Journal**: {pub.journal or 'N/A'} ({pub.year or 'N/A'})")
                lines.append(f"**PMID**: [{pub.pmid}](https://pubmed.ncbi.nlm.nih.gov/{pub.pmid}/)")
                if pub.doi:
                    lines.append(f"**DOI**: {pub.doi}")
                lines.append("")
        else:
            lines.append("No publications found matching this query.")

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)
    except ClientError as exc:
        return f"PubMed search error: {exc}\n\n{DISCLAIMER}"


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True})  # type: ignore[arg-type]
async def get_publication(pmid: str) -> str:
    """Fetch full publication details from PubMed by PMID — title, all authors, abstract, journal, year, DOI, and MeSH terms.

    Use to retrieve the complete metadata and abstract for a specific paper referenced in evidence items or classification guidelines. Essential for verifying the source behind a clinical assertion.

    Use for questions like:
    - "Get the abstract for PMID 27993330 (AMP/ASCO/CAP guidelines)"
    - "What is the paper behind PMID 35101336?"

    Args:
        pmid: PubMed ID (e.g., 27993330, 35101336).
    """
    try:
        pub = await pubmed_client.get_publication(pmid)

        lines = [f"## 📄 Publication — PMID {pub.pmid}\n"]
        lines.append(f"**Title**: {pub.title or 'N/A'}")
        lines.append(f"**Authors**: {', '.join(pub.authors) if pub.authors else 'N/A'}")
        lines.append(f"**Journal**: {pub.journal or 'N/A'} ({pub.year or 'N/A'})")
        lines.append(
            f"**PubMed**: [https://pubmed.ncbi.nlm.nih.gov/{pub.pmid}/]"
            f"(https://pubmed.ncbi.nlm.nih.gov/{pub.pmid}/)"
        )
        if pub.doi:
            lines.append(f"**DOI**: {pub.doi}")
        if pub.abstract:
            lines.append(f"\n### Abstract\n\n{pub.abstract}")
        if pub.mesh_terms:
            lines.append(f"\n**MeSH Terms**: {', '.join(pub.mesh_terms)}")

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)
    except ClientError as exc:
        return f"PubMed fetch error: {exc}\n\n{DISCLAIMER}"


@mcp.tool(annotations={"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True})  # type: ignore[arg-type]
async def predict_variant_effect(
    gene: str,
    variant: str,
    hgvs_id: str | None = None,
) -> str:
    """Aggregate in-silico pathogenicity predictions from 7 computational tools: SIFT, PolyPhen-2, REVEL, CADD, AlphaMissense, GERP++, and phyloP.

    Queries MyVariant.info/dbNSFP and returns individual scores, consensus classification (Damaging/Benign/Mixed), ClinGen SVI-calibrated REVEL strength assessment (Pejaver et al. 2022), and clinical interpretation for ACMG PP3/BP4 and oncogenicity OP1/SBP1 evidence codes. Output includes a structured predictor table ideal for comparison visualizations.

    Use for questions like:
    - "What do computational predictors say about BRAF V600E?"
    - "Is TP53 R175H predicted to be damaging?"
    - "Get REVEL and CADD scores for EGFR L858R"

    Args:
        gene: Gene symbol (e.g., BRAF, TP53, EGFR). Required.
        variant: Variant name (e.g., V600E, R175H, L858R). Required.
        hgvs_id: HGVS genomic ID for direct lookup (e.g., chr7:g.140453136A>T). Optional.
            **Must use GRCh37/hg19 coordinates** — MyVariant.info indexes on hg19.
            Without hgvs_id, searches by gene+protein change (build-agnostic).
    """
    try:
        preds = None
        if hgvs_id:
            preds = await myvariant_client.get_predictions(hgvs_id)
        else:
            preds = await myvariant_client.search_variant(gene, variant)

        if preds is None:
            return (
                f"No in-silico predictions found for {gene} {variant}.\n\n"
                "This may occur for non-coding variants, indels, or novel variants "
                "not yet indexed in dbNSFP.\n\n" + DISCLAIMER
            )

        lines = [f"## 🤖 In-Silico Predictions — {gene} {variant}\n"]
        lines.append(
            f"**Consensus**: {preds.consensus or 'N/A'} "
            f"({preds.damaging_count}/{preds.total_predictors} damaging, "
            f"{preds.benign_count}/{preds.total_predictors} benign)\n"
        )

        lines.append("### Individual Predictor Scores\n")
        lines.append("| Predictor | Score | Prediction | Threshold |")
        lines.append("|-----------|-------|------------|-----------|")

        if preds.sift_score is not None:
            lines.append(
                f"| SIFT | {preds.sift_score:.4f} | {preds.sift_prediction or '—'} | <0.05 = damaging |"
            )
        if preds.polyphen2_score is not None:
            lines.append(
                f"| PolyPhen-2 | {preds.polyphen2_score:.4f} | {preds.polyphen2_prediction or '—'} "
                f"| >0.85 = prob. damaging |"
            )
        if preds.revel_score is not None:
            svi_str = f" **→ {preds.revel_acmg_strength}**" if preds.revel_acmg_strength else ""
            lines.append(
                f"| REVEL | {preds.revel_score:.4f} | —{svi_str} "
                f"| ClinGen SVI: ≥0.773 PP3_strong, ≥0.644 PP3_mod |"
            )
        if preds.cadd_phred is not None:
            lines.append(f"| CADD (phred) | {preds.cadd_phred:.1f} | — | >20 = top 1% |")
        if preds.alphamissense_score is not None:
            lines.append(
                f"| AlphaMissense | {preds.alphamissense_score:.4f} | "
                f"{preds.alphamissense_prediction or '—'} | >0.564 = likely pathogenic |"
            )
        if preds.gerp_score is not None:
            lines.append(f"| GERP++ | {preds.gerp_score:.2f} | — | >2 = conserved |")
        if preds.phylop_score is not None:
            lines.append(f"| phyloP | {preds.phylop_score:.2f} | — | >0 = conserved |")

        lines.append("\n### Clinical Interpretation\n")
        if preds.consensus == "Damaging":
            lines.append("🔴 **Most predictors agree: DAMAGING**")
            lines.append(
                "- Supports **PP3** (ACMG germline: computational evidence supports deleterious)"
            )
            lines.append(
                "- Supports **OP1** (oncogenicity SOP: all computational evidence "
                "supports oncogenic, +1 point)"
            )
        elif preds.consensus == "Benign":
            lines.append("🟢 **Most predictors agree: BENIGN**")
            lines.append(
                "- Supports **BP4** (ACMG germline: computational evidence suggests no impact)"
            )
            lines.append(
                "- Supports **SBP1** (oncogenicity SOP: all computational evidence "
                "suggests benign, −1 point)"
            )
        else:
            lines.append("🟡 **Mixed predictions** — no strong computational consensus")
            lines.append("- PP3/BP4 and OP1/SBP1 should NOT be applied")

        # ClinGen SVI-calibrated REVEL strength (Pejaver et al. 2022)
        if preds.revel_acmg_strength:
            lines.append(
                f"\n**ClinGen SVI REVEL assessment**: {preds.revel_acmg_strength} "
                f"(score {preds.revel_score:.3f}, Pejaver et al. 2022, PMID: 36413997)"
            )

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)
    except ClientError as exc:
        return f"MyVariant.info query error: {exc}\n\n{DISCLAIMER}"


# ============================================================================
# Entry point
# ============================================================================


def main() -> None:
    """Run the MCP server via stdio transport."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
