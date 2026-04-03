"""NCBI E-utilities client for ClinVar."""

from __future__ import annotations

import logging
import os
from typing import Any
from xml.etree import ElementTree

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.constants import (
    CLINVAR_BASE_URL,
    CLINVAR_EFETCH_URL,
    CLINVAR_ESEARCH_URL,
    CLINVAR_ESUMMARY_URL,
    CLINVAR_RATE_LIMIT,
    CLINVAR_REVIEW_STARS,
)
from variant_mcp.models.evidence import ClinVarVariant

logger = logging.getLogger(__name__)


class ClinVarClient(BaseClient):
    """Client for NCBI E-utilities (ClinVar)."""

    def __init__(self) -> None:
        self._api_key = os.environ.get("NCBI_API_KEY")
        rate = 10 if self._api_key else CLINVAR_RATE_LIMIT
        super().__init__(base_url="https://eutils.ncbi.nlm.nih.gov", rate_limit=rate)

    def _base_params(self) -> dict[str, str]:
        params: dict[str, str] = {"db": "clinvar", "retmode": "json"}
        if self._api_key:
            params["api_key"] = self._api_key
        return params

    async def search_variants(
        self,
        *,
        gene: str | None = None,
        variant_name: str | None = None,
        clinical_significance: str | None = None,
        review_status: str | None = None,
        disease: str | None = None,
        limit: int = 25,
    ) -> list[ClinVarVariant]:
        """Search ClinVar and return parsed variant summaries."""
        query_parts: list[str] = []
        if gene:
            query_parts.append(f"{gene}[gene]")
        if variant_name:
            query_parts.append(f"{variant_name}[variant name]")
        if clinical_significance:
            query_parts.append(f"{clinical_significance}[CLNSIG]")
        if review_status:
            query_parts.append(f"{review_status}[RVSTAT]")
        if disease:
            query_parts.append(f"{disease}[Disease/Phenotype]")

        if not query_parts:
            raise ClientError(
                "At least one search parameter is required (gene, variant_name, "
                "clinical_significance, disease)."
            )

        term = " AND ".join(query_parts)
        ids = await self._esearch(term, retmax=limit)
        if not ids:
            return []
        return await self._esummary_variants(ids)

    async def get_variant_detail(self, variation_id: int) -> ClinVarVariant:
        """Get full ClinVar variant record by variation ID using efetch VCV XML.

        NCBI deprecated ``rettype=variation`` in 2019.  The current endpoint
        requires ``rettype=vcv`` together with ``is_variationid=true`` so that
        the numeric variation-ID is accepted directly.
        """
        params = {
            "db": "clinvar",
            "id": str(variation_id),
            "rettype": "vcv",
            "retmode": "xml",
            "is_variationid": "true",
        }
        if self._api_key:
            params["api_key"] = self._api_key

        response = await self._request("GET", CLINVAR_EFETCH_URL, params=params)
        return self._parse_efetch_xml(response.text, variation_id)

    async def get_variant_by_rsid(self, rsid: str) -> list[ClinVarVariant]:
        """Search ClinVar by dbSNP rsID."""
        clean_rsid = rsid.lower().replace("rs", "")
        ids = await self._esearch(f"rs{clean_rsid}[dbsnp id]")
        if not ids:
            return []
        return await self._esummary_variants(ids)

    async def get_variant_by_hgvs(self, hgvs_expression: str) -> list[ClinVarVariant]:
        """Search ClinVar by HGVS expression."""
        ids = await self._esearch(f"{hgvs_expression}[variant name]")
        if not ids:
            return []
        return await self._esummary_variants(ids)

    async def _esearch(self, term: str, retmax: int = 25) -> list[str]:
        """Run esearch and return list of IDs."""
        params = self._base_params()
        params.update({"term": term, "retmax": str(retmax)})
        response = await self._request("GET", CLINVAR_ESEARCH_URL, params=params)
        data = response.json()
        result = data.get("esearchresult", {})
        if "ERROR" in result:
            raise ClientError(f"ClinVar esearch error: {result['ERROR']}")
        return result.get("idlist", [])  # type: ignore[no-any-return]

    async def _esummary_variants(self, ids: list[str]) -> list[ClinVarVariant]:
        """Run esummary for given IDs and parse to ClinVarVariant list."""
        params = self._base_params()
        params["id"] = ",".join(ids)
        response = await self._request("GET", CLINVAR_ESUMMARY_URL, params=params)
        data = response.json()
        result = data.get("result", {})
        uid_list = result.get("uids", [])

        variants = []
        for uid in uid_list:
            doc = result.get(uid, {})
            if not doc:
                continue
            variants.append(self._parse_esummary_doc(doc, uid))
        return variants

    @staticmethod
    def _parse_esummary_doc(doc: dict[str, Any], uid: str) -> ClinVarVariant:
        """Parse an esummary document into a ClinVarVariant.

        ClinVar API (2024+) returns three separate classification fields:
        - germline_classification: germline pathogenicity
        - oncogenicity_classification: somatic oncogenicity
        - clinical_impact_classification: somatic clinical impact (AMP tiers)
        The legacy clinical_significance field is deprecated (returns null).
        """
        genes_list = doc.get("genes", [])
        gene_names = [g.get("symbol", "") for g in genes_list if g.get("symbol")]

        # Try new API fields first, fall back to legacy clinical_significance
        significance = ""
        review_status = ""
        last_eval = ""
        conditions: list[str] = []
        classification_parts: list[str] = []

        for field_name in (
            "germline_classification",
            "oncogenicity_classification",
            "clinical_impact_classification",
        ):
            cls_data = doc.get(field_name)
            if isinstance(cls_data, dict) and cls_data.get("description"):
                desc = cls_data["description"]
                classification_parts.append(desc)
                # Use the first non-empty classification for primary fields
                if not significance:
                    significance = desc
                    review_status = cls_data.get("review_status", "")
                    last_eval = cls_data.get("last_evaluated", "")
                # Collect conditions from each classification's trait_set
                for trait in cls_data.get("trait_set", []):
                    name = trait.get("trait_name", "")
                    if name and name not in conditions:
                        conditions.append(name)

        # Fall back to legacy field if new fields are empty
        if not significance:
            clin_sig = doc.get("clinical_significance", {})
            if isinstance(clin_sig, dict):
                significance = clin_sig.get("description", "")
                review_status = clin_sig.get("review_status", "")
                last_eval = clin_sig.get("last_evaluated", "")
            elif clin_sig:
                significance = str(clin_sig)

        # If multiple classification types exist, combine them
        if len(classification_parts) > 1:
            significance = " | ".join(classification_parts)

        # Fall back to legacy trait_set for conditions if none found
        if not conditions:
            conditions = [
                t.get("trait_name", "") for t in doc.get("trait_set", []) if t.get("trait_name")
            ]

        review_stars = CLINVAR_REVIEW_STARS.get(review_status.lower(), "")

        variation_id = int(uid)
        return ClinVarVariant(
            variation_id=variation_id,
            title=doc.get("title", ""),
            clinical_significance=significance,
            review_status=review_status,
            review_stars=review_stars,
            genes=gene_names,
            variation_type=doc.get("variation_type", ""),
            protein_change=doc.get("protein_change", ""),
            conditions=conditions,
            last_evaluated=last_eval,
            clinvar_url=f"{CLINVAR_BASE_URL}/{variation_id}",
        )

    @staticmethod
    def _parse_efetch_xml(xml_text: str, variation_id: int) -> ClinVarVariant:
        """Parse ClinVar efetch VCV XML response.

        The VCV XML (``rettype=vcv``) uses ``ClassifiedRecord`` (not the
        older ``InterpretedRecord``).  Classifications live under
        ``Classifications/GermlineClassification`` (and siblings for
        somatic / oncogenicity).  HGVS expressions are structured as
        ``NucleotideExpression`` / ``ProteinExpression`` child elements
        with ``sequenceAccessionVersion`` + ``change`` attributes.
        """
        try:
            root = ElementTree.fromstring(xml_text)
        except ElementTree.ParseError as e:
            raise ClientError(f"Failed to parse ClinVar XML: {e}") from e

        variant = ClinVarVariant(
            variation_id=variation_id,
            clinvar_url=f"{CLINVAR_BASE_URL}/{variation_id}",
        )

        # Parse VariationArchive
        va = root.find(".//VariationArchive")
        if va is not None:
            variant.title = va.get("VariationName", "")
            variant.variation_type = va.get("VariationType", "")

        # Parse ClassifiedRecord (VCV format) or fall back to InterpretedRecord
        ir = root.find(".//ClassifiedRecord")
        if ir is None:
            ir = root.find(".//InterpretedRecord")

        if ir is not None:
            # ---- Clinical significance from Classifications block ----
            # VCV XML puts classifications under Classifications/<Type>
            cls_block = ir.find("Classifications")
            classification_parts: list[str] = []
            if cls_block is not None:
                for cls_type in (
                    "GermlineClassification",
                    "SomaticClinicalImpact",
                    "OncogenicityClassification",
                ):
                    cls_el = cls_block.find(cls_type)
                    if cls_el is not None:
                        desc_el = cls_el.find("Description")
                        if desc_el is not None and desc_el.text:
                            classification_parts.append(desc_el.text)
                            # Use the first non-empty classification for primary fields
                            if not variant.clinical_significance:
                                variant.clinical_significance = desc_el.text
                                review_el = cls_el.find("ReviewStatus")
                                if review_el is not None and review_el.text:
                                    variant.review_status = review_el.text
                                    variant.review_stars = CLINVAR_REVIEW_STARS.get(
                                        review_el.text.lower(), ""
                                    )

                if len(classification_parts) > 1:
                    variant.clinical_significance = " | ".join(classification_parts)

            # Fall back to legacy Interpretation/Classification element
            if not variant.clinical_significance:
                interp = ir.find(".//Interpretation") or ir.find(".//Classification")
                if interp is not None:
                    desc = interp.find("Description")
                    if desc is not None and desc.text:
                        variant.clinical_significance = desc.text
                    review = interp.find("ReviewStatus")
                    if review is not None and review.text:
                        variant.review_status = review.text
                        variant.review_stars = CLINVAR_REVIEW_STARS.get(review.text.lower(), "")
                    last_eval_el = interp.find("DateLastEvaluated")
                    if last_eval_el is not None and last_eval_el.text:
                        variant.last_evaluated = last_eval_el.text

            # ---- HGVS expressions ----
            # VCV format: SimpleAllele/HGVSlist/HGVS with child elements
            hgvs_set: list[str] = []
            for hgvs_el in ir.findall(".//SimpleAllele/HGVSlist/HGVS"):
                # Try structured sub-elements first (VCV format)
                for sub_tag in ("NucleotideExpression", "ProteinExpression"):
                    sub = hgvs_el.find(sub_tag)
                    if sub is not None:
                        acc = sub.get("sequenceAccessionVersion", "")
                        change = sub.get("change", "")
                        if acc and change:
                            hgvs_set.append(f"{acc}:{change}")
                # Fall back to text / Expression attribute (legacy)
                if not hgvs_set:
                    text = (hgvs_el.text or "").strip() or hgvs_el.get("Expression", "")
                    if text:
                        hgvs_set.append(text)
            # Also check legacy paths
            if not hgvs_set:
                for hgvs_el in ir.findall(".//HGVSlist/HGVS") + ir.findall(".//HGVSExpression"):
                    text = (hgvs_el.text or "").strip() or hgvs_el.get("Expression", "")
                    if text:
                        hgvs_set.append(text)
            variant.hgvs_expressions = hgvs_set

            # ---- Genes (deduplicated) ----
            gene_names_seen: set[str] = set()
            gene_names: list[str] = []
            for gene_el in ir.findall(".//GeneList/Gene"):
                symbol = gene_el.get("Symbol")
                if symbol and symbol not in gene_names_seen:
                    gene_names_seen.add(symbol)
                    gene_names.append(symbol)
            variant.genes = gene_names

            # ---- Conditions (deduplicated) ----
            conditions_seen: set[str] = set()
            conditions: list[str] = []
            for trait in ir.findall(".//TraitSet/Trait"):
                name_el = trait.find("Name/ElementValue[@Type='Preferred']")
                if name_el is not None and name_el.text:
                    cond = name_el.text
                    if cond not in conditions_seen:
                        conditions_seen.add(cond)
                        conditions.append(cond)
            variant.conditions = conditions

            # ---- Submitter classifications ----
            submitter_list = []
            for scv in ir.findall(".//ClinicalAssertionList/ClinicalAssertion"):
                submitter_el = scv.find(".//ClinVarAccession")
                # VCV: Classification/GermlineClassification (or Description)
                interp_el = scv.find(".//Classification/GermlineClassification")
                if interp_el is None:
                    interp_el = scv.find(".//Classification/Description")
                if interp_el is None:
                    interp_el = scv.find(".//Interpretation/Description")
                if submitter_el is not None:
                    entry = {
                        "submitter": submitter_el.get("SubmitterName", "Unknown"),
                        "accession": submitter_el.get("Accession", ""),
                        "classification": (interp_el.text if interp_el is not None else "N/A"),
                    }
                    submitter_list.append(entry)
            variant.submitter_classifications = submitter_list
            variant.submitter_count = len(submitter_list)

            # Check for conflicting interpretations
            classifications = {
                s["classification"] for s in submitter_list if s["classification"] != "N/A"
            }
            variant.conflicting = len(classifications) > 1

        # rsID
        for xref in root.findall(".//XRef"):
            if xref.get("DB") == "dbSNP":
                variant.rsid = f"rs{xref.get('ID', '')}"
                break

        return variant
