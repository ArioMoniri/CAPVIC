# Claude Code Prompt: Build a Multi-Source Precision Oncology & Variant Classification MCP Server

> **PASTE THIS ENTIRE PROMPT INTO CLAUDE CODE.**
> It builds a production-grade MCP server that integrates 5+ authoritative data sources
> and implements AMP/ASCO/CAP + ClinGen/CGC/VICC classification frameworks.

---

## 🎯 Who This Is For

Researchers building AI models that predict whether genetic variants are **Pathogenic** or **Benign** or so on. This MCP gives any AI assistant (Claude, GPT, etc.) or any other structure the ability to query real clinical evidence, apply formal classification frameworks, and produce structured reports — turning the assistant into a virtual molecular tumor board.

---

## 🏗️ Architecture Overview

**Server name**: `variant_mcp`
**Framework**: FastMCP (Python MCP SDK)
**Transport**: stdio (for Claude Desktop / Claude Code)
**Language**: Python 3.11+
**Data sources** (all free for academic/research use, no auth required except OncoKB):

| Source | API | What It Provides |
|--------|-----|------------------|
| **CIViC** | GraphQL (`civicdb.org/api/graphql`) | Expert-curated clinical evidence: therapeutic, diagnostic, prognostic, oncogenic |
| **ClinVar** | NCBI E-utilities REST (`eutils.ncbi.nlm.nih.gov`) | Aggregate pathogenicity classifications from 2000+ labs worldwide |
| **OncoKB** | REST (`oncokb.org/api/v1`) | FDA-recognized oncogenicity + therapeutic levels (free academic API key) |
| **VICC MetaKB** | REST (`search.cancervariants.org`) | Harmonized interpretations from 6 knowledgebases (CIViC, OncoKB, CGI, JAX-CKB, MMatch, PMKB) |
| **NCBI Entrez / PubMed** | E-utilities | Literature references, gene info, cross-linking |

**Classification frameworks implemented as computational logic (not external APIs)**:

| Framework | Source | Purpose |
|-----------|--------|---------|
| **AMP/ASCO/CAP 4-Tier** (2017, 2025 update) | Li et al., J Mol Diagn 2017 | Somatic variant clinical significance: Tier I (strong) → Tier IV (benign) |
| **ClinGen/CGC/VICC Oncogenicity SOP** | Horak et al., Genet Med 2022 | Point-based oncogenicity: Oncogenic → Benign with evidence codes (OVS1, OS1-3, OM1-4, OP1-4, SBVS1, SBS1-2, SBP1-2) |
| **ACMG/AMP 5-Tier Germline** | Richards et al., Genet Med 2015 | Germline pathogenicity: Pathogenic → Benign with evidence codes (PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7) |

---

## Project Structure

```
variant-mcp/
├── src/
│   └── variant_mcp/
│       ├── __init__.py
│       ├── server.py                  # FastMCP server + all tool registrations
│       ├── clients/
│       │   ├── __init__.py
│       │   ├── civic_client.py        # CIViC V2 GraphQL client
│       │   ├── clinvar_client.py      # NCBI E-utilities client for ClinVar
│       │   ├── oncokb_client.py       # OncoKB REST API client
│       │   ├── metakb_client.py       # VICC MetaKB search client
│       │   └── base_client.py         # Shared async HTTP + rate limiting
│       ├── classification/
│       │   ├── __init__.py
│       │   ├── amp_asco_cap.py        # AMP/ASCO/CAP 4-tier somatic classifier
│       │   ├── oncogenicity_sop.py    # ClinGen/CGC/VICC oncogenicity SOP scorer
│       │   └── acmg_amp.py            # ACMG/AMP germline pathogenicity helper
│       ├── models/
│       │   ├── __init__.py
│       │   ├── inputs.py             # All Pydantic input models
│       │   ├── evidence.py           # Unified evidence data models
│       │   └── classification.py     # Classification result models
│       ├── formatters/
│       │   ├── __init__.py
│       │   ├── reports.py            # Markdown report generators
│       │   └── tables.py             # Summary table formatters
│       ├── queries/
│       │   ├── __init__.py
│       │   └── civic_graphql.py      # CIViC GraphQL query strings
│       └── constants.py              # URLs, enums, disclaimers, evidence codes
├── tests/
│   ├── test_clients.py
│   ├── test_classification.py
│   ├── test_formatters.py
│   └── test_models.py
├── Dockerfile
├── docker-compose.yml
├── pyproject.toml
├── README.md
└── .env.example
```

Dependencies in `pyproject.toml`:
```toml
[project]
name = "variant-mcp"
version = "1.0.0"
requires-python = ">=3.11"
dependencies = [
    "mcp[cli]>=1.0.0",
    "httpx>=0.27.0",
    "pydantic>=2.0",
    "lxml>=5.0",          # For parsing ClinVar XML responses
]
```

---

## PHASE 1: Constants & Data Models

### `constants.py`

```python
# API Endpoints
CIVIC_GRAPHQL_URL = "https://civicdb.org/api/graphql"
CIVIC_BASE_URL = "https://civicdb.org"
CLINVAR_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CLINVAR_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CLINVAR_BASE_URL = "https://www.ncbi.nlm.nih.gov/clinvar/variation"
ONCOKB_API_URL = "https://www.oncokb.org/api/v1"
METAKB_SEARCH_URL = "https://search.cancervariants.org"  # check actual API endpoint

# Rate limits
CIVIC_RATE_LIMIT = 3      # 3 req/sec without API key
CLINVAR_RATE_LIMIT = 3     # 3 req/sec without API key, 10 with
ONCOKB_RATE_LIMIT = 5      # varies by license

# Timeouts
REQUEST_TIMEOUT = 30

# === AMP/ASCO/CAP 4-Tier Classification (2017 + 2025 update) ===
AMP_TIER_DEFINITIONS = {
    "I": {
        "name": "Tier I — Strong Clinical Significance",
        "level_A": "FDA-approved therapy, included in professional guidelines (NCCN, CAP, ASCO, ESMO), or well-powered studies with consensus",
        "level_B": "Well-powered studies with expert consensus from clinical trials or large cohorts",
    },
    "II": {
        "name": "Tier II — Potential Clinical Significance",
        "level_C": "FDA-approved therapies for different tumor types, or investigational therapies with evidence from clinical trials",
        "level_D": "Preclinical trials or case reports showing potential association with diagnosis, prognosis, or therapy",
    },
    "III": {
        "name": "Tier III — Unknown Clinical Significance",
        "description": "Variant not observed at a significant frequency in population databases, no convincing evidence of cancer association published",
    },
    "IV": {
        "name": "Tier IV — Benign or Likely Benign",
        "description": "Observed at significant frequency in population databases (gnomAD/ExAC). No published evidence of cancer association.",
    }
}

# === ClinGen/CGC/VICC Oncogenicity SOP Evidence Codes ===
# Reference: Horak et al., Genetics in Medicine 2022
ONCOGENICITY_EVIDENCE_CODES = {
    # Oncogenic codes (positive evidence)
    "OVS1": {"strength": "very_strong", "points": 8, "description": "Null variant in bona fide tumor suppressor gene"},
    "OS1":  {"strength": "strong", "points": 4, "description": "Same amino acid change as established oncogenic variant"},
    "OS2":  {"strength": "strong", "points": 4, "description": "Well-established functional studies show oncogenic effect"},
    "OS3":  {"strength": "strong", "points": 4, "description": "Located in hotspot with sufficient statistical evidence"},
    "OM1":  {"strength": "moderate", "points": 2, "description": "Located in critical/functional domain without benign variation"},
    "OM2":  {"strength": "moderate", "points": 2, "description": "Protein length changes in known oncogene/TSG"},
    "OM3":  {"strength": "moderate", "points": 2, "description": "At same position as established oncogenic missense"},
    "OM4":  {"strength": "moderate", "points": 2, "description": "Absent/extremely rare in population databases"},
    "OP1":  {"strength": "supporting", "points": 1, "description": "All computational evidence supports oncogenic effect"},
    "OP2":  {"strength": "supporting", "points": 1, "description": "Somatic variant in gene with single cancer etiology"},
    "OP3":  {"strength": "supporting", "points": 1, "description": "Located in hotspot with limited statistical evidence"},
    "OP4":  {"strength": "supporting", "points": 1, "description": "Absent in population databases"},
    # Benign codes (negative evidence)
    "SBVS1": {"strength": "very_strong_benign", "points": -8, "description": "MAF >5% in any general continental population in gnomAD"},
    "SBS1":  {"strength": "strong_benign", "points": -4, "description": "MAF >1% in any general continental population in gnomAD"},
    "SBS2":  {"strength": "strong_benign", "points": -4, "description": "Well-established functional studies show no oncogenic effect"},
    "SBP1":  {"strength": "supporting_benign", "points": -1, "description": "All computational evidence suggests benign"},
    "SBP2":  {"strength": "supporting_benign", "points": -1, "description": "Synonymous change with no predicted splice impact"},
}

# Oncogenicity classification thresholds (point-based)
ONCOGENICITY_THRESHOLDS = {
    "Oncogenic": 10,          # >= 10 points
    "Likely Oncogenic": 6,     # >= 6 and < 10 points
    "VUS": None,               # between -1 and 5 (inclusive)
    "Likely Benign": -6,       # <= -6 and > -10
    "Benign": -10,             # <= -10
}

# === ACMG/AMP Germline Classification Terms ===
ACMG_CLASSIFICATIONS = [
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
]

# ClinVar review status star ratings
CLINVAR_REVIEW_STARS = {
    "practice guideline": "★★★★",
    "reviewed by expert panel": "★★★",
    "criteria provided, multiple submitters, no conflicts": "★★",
    "criteria provided, single submitter": "★",
    "no assertion criteria provided": "☆",
    "no assertion provided": "☆",
}

# CIViC Evidence Level descriptions
CIVIC_EVIDENCE_LEVELS = {
    "A": "Validated — Proven/consensus in human medicine",
    "B": "Clinical — Clinical trial or primary patient data",
    "C": "Case Study — Individual case reports",
    "D": "Preclinical — In vivo/in vitro models",
    "E": "Inferential — Indirect/computational",
}

# OncoKB Therapeutic Levels of Evidence
ONCOKB_LEVELS = {
    "LEVEL_1": "FDA-recognized biomarker predictive of response in this indication",
    "LEVEL_2": "Standard care biomarker (NCCN/expert panels) in this indication",
    "LEVEL_3A": "Compelling clinical evidence supports biomarker as predictive in this indication",
    "LEVEL_3B": "Standard care or investigational biomarker predictive in another indication",
    "LEVEL_4": "Compelling biological evidence supports biomarker as predictive",
    "LEVEL_R1": "Standard care biomarker predictive of resistance (FDA/guidelines)",
    "LEVEL_R2": "Compelling clinical evidence supports biomarker as predictive of resistance",
}

DISCLAIMER = """⚠️ **RESEARCH DISCLAIMER**: This tool aggregates data from CIViC, ClinVar, OncoKB, and VICC MetaKB for **research and educational purposes only**. Classification suggestions follow published AMP/ASCO/CAP, ClinGen/CGC/VICC, and ACMG/AMP frameworks but are **computationally derived approximations** — NOT validated clinical interpretations. Clinical variant classification MUST be performed by qualified molecular pathologists using validated laboratory procedures. Do not use for clinical decision-making without professional review."""
```

### `models/inputs.py`

Define Pydantic v2 models for ALL tools. Each model MUST have:
- `model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')`
- Descriptive Field() with min/max constraints
- Enum types for controlled vocabularies

Key input models to create:
- `SearchEvidenceInput` — multi-source evidence search (disease, gene, variant, therapy, evidence_type, source_databases[], limit)
- `ClassifyVariantInput` — full classification (gene, variant, disease/cancer_type, variant_origin [somatic/germline])
- `GetClinVarInput` — ClinVar lookup (gene, variant, variation_id, clinical_significance)
- `SearchCIViCInput` — CIViC-specific search
- `AnnotateOncoKBInput` — OncoKB annotation (gene, variant, tumor_type)
- `SearchMetaKBInput` — VICC MetaKB harmonized search
- `LookupGeneInput` / `LookupDiseaseInput` / `LookupTherapyInput` — typeahead inputs
- `AMPTierInput` — AMP/ASCO/CAP tier assignment
- `OncogenicityScoringInput` — ClinGen/CGC/VICC evidence code scoring

---

## PHASE 2: API Clients

### `clients/base_client.py`
Shared infrastructure:
- Async `httpx.AsyncClient` with configurable timeout
- Rate-limiter using `asyncio.Semaphore` + `asyncio.sleep`
- Retry logic with exponential backoff (max 3 retries)
- Unified error handling returning actionable messages
- JSON + XML response parsing helpers

### `clients/civic_client.py` — CIViC V2 GraphQL
Endpoint: `POST https://civicdb.org/api/graphql`
No auth required. Rate limit: 3 req/sec.

Methods:
- `search_evidence(disease, gene, variant, therapy, evidence_type, significance, first, after)` — returns evidence items with pagination
- `get_gene(name)` — gene details + variants
- `get_variant(variant_id)` — variant + molecular profiles
- `get_evidence_item(evidence_id)` — full detail
- `search_assertions(disease, gene, therapy, significance, first, after)` — curated assertions
- `search_typeahead(entity_type, query_term)` — autocomplete for genes/diseases/therapies

Use GraphQL queries from `queries/civic_graphql.py`. Key fields to always request on evidence items:
`id, name, status, evidenceType, evidenceLevel, evidenceDirection, evidenceRating, significance, description, therapyInteractionType, molecularProfile { id name variants { name id } }, disease { id name doid }, therapies { id name ncitId }, phenotypes { id name hpoId }, source { id citation sourceUrl citationId sourceType }`

### `clients/clinvar_client.py` — NCBI E-utilities
Endpoints: esearch, esummary, efetch at `eutils.ncbi.nlm.nih.gov`
No auth required (3 req/sec). Optional NCBI API key for 10 req/sec (env var `NCBI_API_KEY`).

Methods:
- `search_variants(gene, variant_name, clinical_significance, review_status, disease)` — esearch then esummary
  - Build query like: `{gene}[gene] AND {significance}[CLNSIG]`
  - Parse esummary JSON for: variation_id, title, clinical_significance, review_status, genes, variation_type, protein_change, last_evaluated
- `get_variant_detail(variation_id)` — efetch with `rettype=variation` XML
  - Parse XML to extract: all submitter classifications, review status, conditions, HGVS expressions, dbSNP rsID, genomic coordinates, submitter details with evidence
- `get_variant_by_rsid(rsid)` — esearch with rsID then detail
- `get_variant_by_hgvs(hgvs_expression)` — esearch with HGVS then detail

Key ClinVar search fields: `[gene]`, `[CLNSIG]` (clinical significance), `[RVSTAT]` (review status), `[Disease/Phenotype]`, `[VRTP]` (variant type), `[chr]`, `[chrpos37]`/`[chrpos38]`

### `clients/oncokb_client.py` — OncoKB REST API
Endpoint: `https://www.oncokb.org/api/v1`
Requires API token (free for academic). Env var: `ONCOKB_API_TOKEN`. If not set, tool should gracefully degrade and explain how to get a free academic token.

Methods:
- `annotate_mutation(gene, variant, tumor_type)` — `/annotate/mutations/byProteinChange`
  - Returns: oncogenic effect, mutation effect, therapeutic levels, diagnostic/prognostic implications
- `annotate_by_hgvsg(hgvsg, tumor_type, ref_genome)` — `/annotate/mutations/byHGVSg`
- `get_gene_info(gene)` — gene summary, oncogene/TSG status
- `get_cancer_gene_list()` — all OncoKB cancer genes
- `get_all_actionable_variants()` — therapeutic biomarkers

Key OncoKB response fields: `oncogenic` (Oncogenic, Likely Oncogenic, Predicted Oncogenic, Likely Neutral, Inconclusive, Unknown), `mutationEffect.knownEffect` (Gain-of-function, Loss-of-function, etc.), `highestSensitiveLevel`, `highestResistanceLevel`, `treatments[]`, `geneSummary`, `variantSummary`

### `clients/metakb_client.py` — VICC MetaKB
Endpoint: `https://search.cancervariants.org` (check for REST API or GraphQL endpoint)
No auth required.

Methods:
- `search(gene, variant, disease)` — harmonized search across all member knowledgebases
  - Returns interpretations from CIViC, OncoKB, CGI, JAX-CKB, MMatch, PMKB
  - Each result includes: source, disease, drugs, evidence_level, clinical_significance

If the MetaKB REST API is not directly accessible at the time of implementation, implement a fallback that queries CIViC + ClinVar and notes that MetaKB web search is available at search.cancervariants.org for manual cross-reference.

---

## PHASE 3: Classification Engines

### `classification/amp_asco_cap.py` — AMP/ASCO/CAP 4-Tier Somatic Classifier

Implement the 2017 guidelines with 2025 update awareness:

```python
class AMPTierClassifier:
    """
    AMP/ASCO/CAP somatic variant clinical significance tiering.

    Tier I: Strong clinical significance
      Level A: FDA-approved therapy, guideline-recommended
      Level B: Well-powered clinical studies with consensus
    Tier II: Potential clinical significance
      Level C: FDA-approved for different tumor type, or investigational in clinical trials
      Level D: Preclinical/case report evidence
    Tier III: Unknown clinical significance
    Tier IV: Benign or Likely Benign

    Input: aggregated evidence from CIViC, ClinVar, OncoKB
    Output: suggested tier with evidence summary and confidence
    """

    def classify(self, evidence_bundle: EvidenceBundle) -> AMPTierResult:
        """
        Logic:
        1. Check OncoKB levels: L1 → Tier I/A, L2 → Tier I/B, L3A → Tier II/C, L3B/L4 → Tier II/D
        2. Check CIViC evidence: Level A accepted → Tier I, Level B → Tier I/B or Tier II, etc.
        3. Check CIViC assertions (higher confidence than evidence items)
        4. Cross-reference ClinVar: Benign/Likely Benign → supports Tier IV
        5. Population frequency (if available): >1% gnomAD → Tier IV
        6. If no evidence found → Tier III
        7. Return tier with evidence trail showing which sources contributed
        """
```

### `classification/oncogenicity_sop.py` — ClinGen/CGC/VICC Oncogenicity Scorer

Implement the point-based system from Horak et al., Genet Med 2022:

```python
class OncogenicityScorer:
    """
    ClinGen/CGC/VICC Standard Operating Procedure for classification
    of oncogenicity of somatic variants.

    Uses a point-based system with evidence codes:
    Oncogenic: OVS1, OS1-3, OM1-4, OP1-4
    Benign: SBVS1, SBS1-2, SBP1-2

    Classification thresholds:
    >= 10 points → Oncogenic
    6-9 points  → Likely Oncogenic
    0-5 points  → VUS
    -1 to -5    → VUS
    -6 to -9    → Likely Benign
    <= -10      → Benign
    """

    def score_variant(self, gene: str, variant: str, evidence_codes: list[str]) -> OncogenicityScoringResult:
        """
        Accept explicit evidence codes OR attempt to auto-detect applicable codes
        from available data:
        - OVS1: check if null variant (nonsense, frameshift, splice) in known TSG
        - OS3: check CIViC/OncoKB hotspot annotations
        - OM4/OP4: check population frequency from ClinVar gnomAD data
        - OP1: check in silico predictions if available
        - SBVS1/SBS1: check gnomAD frequency > 5% or > 1%

        Returns: total points, classification, evidence code breakdown, confidence
        """

    def explain_codes(self, codes: list[str]) -> str:
        """Return human-readable explanation of which evidence codes apply and why."""
```

### `classification/acmg_amp.py` — ACMG/AMP Germline Pathogenicity Helper

Lighter implementation — provides reference information and helps interpret ClinVar data:

```python
class ACMGAMPHelper:
    """
    Reference helper for ACMG/AMP germline variant pathogenicity guidelines.
    Richards et al., Genetics in Medicine 2015.

    NOTE: This does NOT perform full ACMG classification (which requires
    manual expert review). It helps teams understand ClinVar classifications
    and the evidence criteria framework.

    Evidence criteria groups:
    - Very Strong: PVS1 (null variant in gene with known LOF mechanism)
    - Strong: PS1-4 (same AA, functional, segregation, prevalence)
    - Moderate: PM1-6 (functional domain, rare, in-trans, protein impact, cosegregation, variant type)
    - Supporting: PP1-5 (cosegregation, functional, computational, phenotype, reputable source)
    - Benign standalone: BA1 (allele freq > 5%)
    - Benign strong: BS1-4 (freq, functional, segregation, reputable source)
    - Benign supporting: BP1-7 (missense in truncation gene, functional, in-trans, silent, etc.)
    """

    def explain_clinvar_classification(self, clinvar_data: dict) -> str:
        """
        Given ClinVar esummary data, explain:
        - What the aggregate classification means
        - What review status stars mean
        - Whether there are conflicting interpretations
        - Which submitters (labs/expert panels) contributed
        - Date of last evaluation
        """

    def get_criteria_reference(self, criteria_code: str) -> str:
        """Return description and typical evidence for a given ACMG code like PVS1, PM2, etc."""
```

---

## PHASE 4: MCP Tools (server.py)

Register ALL tools on `mcp = FastMCP("variant_mcp", instructions="...")`.

### CATEGORY 1: Multi-Source Search Tools (the power tools)

#### Tool 1: `variant_search_evidence`
**THE PRIMARY TOOL.** Queries ALL available data sources simultaneously and returns a unified report.
- Input: gene (required), variant (optional), disease (optional), therapy (optional), evidence_type (optional), sources (optional, default=all)
- Logic:
  1. Query CIViC for clinical evidence
  2. Query ClinVar for pathogenicity classifications
  3. Query OncoKB for oncogenicity annotation (if token available)
  4. Aggregate and deduplicate results
  5. Format unified report
- Output: Multi-source evidence report with summary statistics, cross-database concordance, top evidence items, citations, disclaimer
- Annotations: `readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=True`

#### Tool 2: `variant_classify`
**THE CLASSIFICATION TOOL.** Given a variant, applies appropriate classification framework(s) and returns a structured assessment.
- Input: gene (required), variant (required), disease/cancer_type (optional), variant_origin (required: "somatic" or "germline")
- Logic:
  - If somatic:
    1. Gather evidence from CIViC + ClinVar + OncoKB
    2. Apply AMP/ASCO/CAP 4-tier classification
    3. Apply ClinGen/CGC/VICC Oncogenicity SOP scoring
    4. Return both classifications side-by-side
  - If germline:
    1. Query ClinVar for classifications + review status
    2. Provide ACMG/AMP framework reference
    3. Summarize submitter agreement/conflicts
    4. Return structured interpretation guide
- Output: Dual-framework classification report with evidence trail
- **This is the competition-critical tool** — directly helps teams validate pathogenic vs benign predictions

#### Tool 3: `variant_compare_sources`
Cross-reference the same variant across all databases and highlight concordance/discordance.
- Input: gene, variant
- Output: Side-by-side comparison showing what each source says about the variant

### CATEGORY 2: CIViC-Specific Tools

#### Tool 4: `civic_search_evidence`
Direct CIViC evidence search with all filters (disease, gene, variant, therapy, evidence_type, significance, phenotype, limit).

#### Tool 5: `civic_search_assertions`
Search CIViC curated assertions (higher confidence than raw evidence items).

#### Tool 6: `civic_get_gene` / Tool 7: `civic_get_variant` / Tool 8: `civic_get_evidence_item`
Detail lookups by ID.

### CATEGORY 3: ClinVar Tools

#### Tool 9: `clinvar_search`
Search ClinVar by gene + variant name, clinical significance, review status, disease.
- Returns: variant classifications, review stars, submitter counts, conflicting interpretations flag, HGVS expressions, last evaluation date

#### Tool 10: `clinvar_get_variant`
Get full ClinVar variant record by variation ID, rsID, or HGVS.
- Returns: comprehensive record including ALL submitter classifications, conditions, review statuses

### CATEGORY 4: OncoKB Tools

#### Tool 11: `oncokb_annotate`
Annotate a somatic mutation with OncoKB (requires API token).
- Input: gene, variant (protein change), tumor_type (optional, OncoTree code)
- Returns: oncogenicity, mutation effect, therapeutic levels, treatment options, gene summary
- If no token: return helpful message explaining how to register for free academic access at oncokb.org

#### Tool 12: `oncokb_cancer_genes`
Get the OncoKB cancer gene list (oncogene/TSG classification).
- Useful for teams to check if a gene is a known oncogene or tumor suppressor

### CATEGORY 5: Classification Framework Tools

#### Tool 13: `classify_amp_tier`
Apply AMP/ASCO/CAP tiering to a somatic variant.
- Input: gene, variant, disease, evidence_bundle (or auto-fetch)
- Output: Suggested tier (I-IV) with evidence level (A-D), supporting evidence trail, confidence assessment

#### Tool 14: `score_oncogenicity`
Apply ClinGen/CGC/VICC oncogenicity SOP.
- Input: gene, variant, evidence_codes (optional — if not provided, attempt auto-detection from CIViC/ClinVar data)
- Output: Point score, classification (Oncogenic/Likely Oncogenic/VUS/Likely Benign/Benign), applied evidence codes with explanations

#### Tool 15: `explain_acmg_criteria`
Reference tool for ACMG/AMP germline classification criteria.
- Input: criteria code (e.g., "PVS1", "PM2") OR general query
- Output: Explanation of criteria, typical evidence, examples

### CATEGORY 6: Discovery Tools

#### Tool 16: `lookup_gene` — Gene name autocomplete across all sources
#### Tool 17: `lookup_disease` — Disease name autocomplete (DOID, OncoTree)
#### Tool 18: `lookup_therapy` — Therapy/drug name autocomplete

### CATEGORY 7: Utility Tools

#### Tool 19: `get_classification_frameworks_reference`
Returns a comprehensive reference document explaining:
- AMP/ASCO/CAP 4-tier system with all evidence levels
- ClinGen/CGC/VICC oncogenicity SOP with all evidence codes
- ACMG/AMP 5-tier germline system
- How these frameworks relate to each other
- When to use which framework (somatic vs germline context)

This tool helps teams understand the classification landscape without searching external docs.

#### Tool 20: `variant_pathogenicity_summary`
**Competition-specific.** Given gene + variant, produce a structured pathogenic vs benign evidence summary:
```
## Pathogenic/Oncogenic Evidence
- ClinVar: 15 submitters say Pathogenic, 2 say Likely Pathogenic (★★★)
- CIViC: 8 evidence items support oncogenic function (3 Level A, 5 Level B)
- OncoKB: Oncogenic, Gain-of-function
- Oncogenicity SOP: 12 points → Oncogenic (OVS1 + OS3 + OM4)

## Benign Evidence
- ClinVar: 0 submitters say Benign
- CIViC: 1 evidence item with Benign significance (Level D)
- gnomAD: Not found (supports OM4/OP4)

## Consensus Assessment
STRONG PATHOGENIC — 17/18 ClinVar submitters agree, OncoKB confirms,
CIViC evidence overwhelmingly supports, Oncogenicity SOP score 12/10 threshold.
Confidence: HIGH
```

---

## PHASE 5: Response Formatting

### Every response from every tool MUST include:

1. **Source attribution**: Which database(s) contributed to the response
2. **Evidence quality indicators**: ClinVar stars, CIViC evidence levels/ratings, OncoKB levels
3. **Direct links**: URLs to CIViC evidence, ClinVar variation pages, OncoKB entries
4. **Citations**: PMIDs with author/year when available
5. **Framework references**: When classification is applied, cite the specific guideline paper
6. **Timestamp**: When the data was retrieved (ClinVar updates weekly, CIViC continuously)
7. **DISCLAIMER**: The full research disclaimer on EVERY response

### Report format example for `variant_classify`:

```markdown
# 🧬 Variant Classification Report

## Variant: BRAF V600E
**Gene**: BRAF | **Variant**: p.Val600Glu | **Origin**: Somatic | **Disease Context**: Melanoma

---

## 📊 Multi-Source Evidence Summary

### CIViC (civicdb.org)
- **Total evidence items**: 47 (42 accepted, 5 submitted)
- **Evidence types**: Predictive (35), Prognostic (8), Diagnostic (4)
- **Key finding**: Strong evidence for sensitivity to vemurafenib, dabrafenib, encorafenib

### ClinVar (ncbi.nlm.nih.gov/clinvar)
- **Variation ID**: 13961
- **Germline classification**: Pathogenic/Likely pathogenic (★★★ — expert panel reviewed)
- **Somatic clinical impact**: Tier I (Strong)
- **Oncogenicity**: Oncogenic
- **Submitters**: 45 (42 Pathogenic, 3 Likely Pathogenic — NO conflicts)

### OncoKB (oncokb.org)
- **Oncogenic**: Oncogenic
- **Mutation Effect**: Gain-of-function
- **Highest Sensitive Level**: Level 1 (Melanoma: vemurafenib, dabrafenib, encorafenib ± binimetinib/cobimetinib/trametinib)
- **Resistance**: Level R1 (Resistance to some therapies in colorectal cancer)

---

## 🏷️ AMP/ASCO/CAP Classification

**Tier I — Strong Clinical Significance (Level A)**

Evidence trail:
- FDA-approved targeted therapy (vemurafenib, Zelboraf®) for BRAF V600E melanoma
- NCCN Guidelines recommend BRAF inhibitor ± MEK inhibitor
- Multiple Phase III trials (BRIM-3, coBRIM, COLUMBUS) confirm response

---

## 🔬 ClinGen/CGC/VICC Oncogenicity Assessment

**Classification: ONCOGENIC (Score: 14 points)**

| Code | Points | Evidence |
|------|--------|----------|
| OS1  | +4     | Same amino acid change as established oncogenic variant (V600E is the canonical BRAF activating mutation) |
| OS2  | +4     | Extensive functional studies show constitutive kinase activation |
| OS3  | +4     | Located in BRAF V600 hotspot (>90% of BRAF mutations in melanoma) |
| OM4  | +2     | Extremely rare in gnomAD (absent from controls) |
| **Total** | **14** | **Oncogenic (threshold: ≥10)** |

---

## Citations
1. Li et al. (2017). "Standards and Guidelines for the Interpretation and Reporting..." J Mol Diagn. PMID: 27993330
2. Horak et al. (2022). "Standards for the classification of pathogenicity of somatic variants..." Genet Med. PMID: 35101336
3. Chapman et al. (2011). "Improved survival with vemurafenib in melanoma..." NEJM. PMID: 21639808

---
⚠️ **RESEARCH DISCLAIMER**: [full disclaimer text]
```

---

## PHASE 6: Server Entry Point

```python
mcp = FastMCP(
    "variant_mcp",
    instructions="""Multi-source precision oncology and variant classification MCP server.

    Integrates CIViC, ClinVar, OncoKB, and VICC MetaKB to provide comprehensive
    clinical evidence for genetic variant interpretation. Implements AMP/ASCO/CAP
    somatic tiering, ClinGen/CGC/VICC oncogenicity scoring, and ACMG/AMP germline
    pathogenicity reference.

    KEY TOOLS:
    - variant_search_evidence: Multi-source evidence search (start here)
    - variant_classify: Apply classification frameworks (somatic or germline)
    - variant_pathogenicity_summary: Pathogenic vs Benign evidence comparison
    - score_oncogenicity: ClinGen/CGC/VICC point-based scoring
    - classify_amp_tier: AMP/ASCO/CAP somatic tiering
    - civic_search_evidence / clinvar_search / oncokb_annotate: Source-specific queries
    - lookup_gene / lookup_disease / lookup_therapy: Discovery/autocomplete
    - get_classification_frameworks_reference: Understand the guidelines

    IMPORTANT: All outputs are for research/education only, not clinical decisions.
    OncoKB features require a free academic API token (set ONCOKB_API_TOKEN env var)."""
)
```

---

## PHASE 7: Dockerfile + docker-compose.yml

```dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY pyproject.toml .
RUN pip install --no-cache-dir .
COPY src/ src/
ENV PYTHONPATH=/app/src
CMD ["python", "-m", "variant_mcp.server"]
```

```yaml
version: "3.8"
services:
  variant-mcp:
    build: .
    stdin_open: true
    environment:
      - ONCOKB_API_TOKEN=${ONCOKB_API_TOKEN:-}
      - NCBI_API_KEY=${NCBI_API_KEY:-}
```

`.env.example`:
```
# Optional: OncoKB API token (free for academic use at https://www.oncokb.org/apiAccess)
ONCOKB_API_TOKEN=

# Optional: NCBI API key (increases rate limit from 3 to 10 req/sec)
# Get one at https://www.ncbi.nlm.nih.gov/account/settings/
NCBI_API_KEY=
```

---

## PHASE 8: README.md

Write comprehensive README with:

1. **Project overview** — what it does, which databases, which classification frameworks
2. **Architecture diagram** (ASCII art showing: User ↔ Claude ↔ MCP Server ↔ [CIViC | ClinVar | OncoKB | MetaKB])
3. **Quick start** — `pip install -e .` → `python -m variant_mcp.server`
4. **Docker** — `docker compose build && docker compose up`
5. **Claude Desktop config** — JSON snippet for `claude_desktop_config.json`
6. **Environment variables** — ONCOKB_API_TOKEN, NCBI_API_KEY
7. **Tool reference table** — all 20 tools with descriptions and example inputs
8. **Classification frameworks** — concise explanation of AMP/ASCO/CAP, Oncogenicity SOP, ACMG/AMP
9. **Example workflows**:
   - "Is BRAF V600E pathogenic in melanoma?" → `variant_classify(gene="BRAF", variant="V600E", disease="Melanoma", variant_origin="somatic")`
   - "What therapies work for KRAS G12C in NSCLC?" → `variant_search_evidence(gene="KRAS", variant="G12C", disease="Non-small Cell Lung Cancer")`
   - "Is BRCA1 R1699W pathogenic?" → `variant_classify(gene="BRCA1", variant="R1699W", variant_origin="germline")`
   - "Score TP53 R175H for oncogenicity" → `score_oncogenicity(gene="TP53", variant="R175H")`
   - "Compare what all databases say about PIK3CA H1047R" → `variant_compare_sources(gene="PIK3CA", variant="H1047R")`
10. **Competition use case** — how teams can use this to validate their ML models
11. **Data freshness** — CIViC (live), ClinVar (weekly), OncoKB (continuously updated)
12. **Key references**:
    - Li et al. (2017) AMP/ASCO/CAP Guidelines. PMID: 27993330
    - Horak et al. (2022) ClinGen/CGC/VICC Oncogenicity SOP. PMID: 35101336
    - Richards et al. (2015) ACMG/AMP Germline Guidelines. PMID: 25741868
    - Griffith et al. (2017) CIViC. PMID: 28138153
    - Wagner et al. (2020) VICC Meta-Knowledgebase. PMID: 32246132
13. **Disclaimer**

---

## PHASE 9: Tests

Using `pytest` + `pytest-asyncio`:

1. **`test_models.py`**: Validate all Pydantic models accept/reject correctly
2. **`test_classification.py`**: 
   - Test AMP tier classifier with known variants (BRAF V600E → Tier I, novel VUS → Tier III)
   - Test oncogenicity scorer with example variants from Horak et al. paper
   - Test ACMG helper explains criteria correctly
3. **`test_formatters.py`**: Test report formatters produce expected Markdown
4. **`test_clients.py`**: Mock httpx responses for each client

---

## CRITICAL REQUIREMENTS

1. **DISCLAIMER on every tool response** — no exceptions, no shortcuts
2. **All network calls async** with proper timeout + retry
3. **Graceful degradation**: if OncoKB token not set, tools still work but note the limitation. If a database is temporarily unavailable, report partial results from other sources.
4. **No medical advice framing** — always "evidence suggests" / "classification frameworks indicate", never "the variant IS pathogenic"
5. **Evidence trail transparency** — every classification must show which sources and which evidence codes contributed. No black-box outputs.
6. **ClinVar review status always shown** — this is critical for reliability assessment (★★★★ expert panel >> ☆ no criteria)
7. **OncoKB licensing compliance** — clearly note it's free for academic research only, mention licensing for commercial use
8. **stderr logging only** — we use stdio transport
9. **Error messages must be actionable** — not just "error occurred" but "CIViC returned 0 results for gene 'KRAS'. Try searching without the variant filter, or check the gene symbol with lookup_gene."
10. **Include guideline paper PMIDs** whenever a classification framework is applied

---

## BUILD ORDER

Execute phases in this exact sequence:
1. `constants.py` → `models/` → `queries/civic_graphql.py`
2. `clients/base_client.py` → all 4 clients
3. `classification/` (all 3 modules)
4. `formatters/`
5. `server.py` (register all 20 tools)
6. `Dockerfile` + `docker-compose.yml` + `.env.example`
7. `pyproject.toml`
8. `README.md`
9. Tests
10. Verify: `python -c "from variant_mcp.server import mcp; print('Server OK')"` and `python -m py_compile src/variant_mcp/server.py`

**After each phase, verify compilation. After all phases, run tests.**

---

## PHASE 10: Scientific Enhancements (Bioinformatics Deep Dive)

*Added by a bioinformatics/genetics review of the codebase. These features address critical gaps that a molecular pathologist or genetic counselor would immediately notice.*

### 10a. gnomAD Population Frequency Client 🧬

**Why**: The oncogenicity SOP codes SBVS1 (MAF >5%), SBS1 (MAF >1%), OM4 (absent/rare), and OP4 (absent) all require population allele frequency data. ACMG PM2 (absent from controls) and BA1 (>5%) also depend on this. Currently the system guesses from ClinVar benign labels — a real clinical interpretation requires actual gnomAD frequencies.

**Implementation**: `clients/gnomad_client.py`
- Query gnomAD GraphQL API (`https://gnomad.broadinstitute.org/api`)
- Fetch per-population allele frequencies (global, AFR, AMR, ASJ, EAS, FIN, NFE, SAS)
- Return: allele count, allele number, allele frequency, homozygote count, filtering status
- Parse both gnomAD v4 (GRCh38) and v2 (GRCh37) for backward compatibility

**New model**: `GnomADFrequency` in `models/evidence.py`
```python
class GnomADFrequency(BaseModel):
    allele_frequency: float | None = None       # Global AF
    allele_count: int | None = None
    allele_number: int | None = None
    homozygote_count: int | None = None
    population_frequencies: dict[str, float] = Field(default_factory=dict)  # Per-population AFs
    genome_version: str = "GRCh38"
    filter_status: str | None = None             # PASS, AC0, InbreedingCoeff, etc.
    source: str = "gnomAD v4"
```

**Integration**: Update `OncogenicityScorer._auto_detect_codes()` and `ACMGAMPHelper` to use actual frequencies.

### 10b. HGVS Notation Parser & Variant Normalizer 🔀

**Why**: Clinicians type "BRAF V600E", databases need "NM_004333.6:c.1799T>A" (ClinVar) or "p.Val600Glu" (HGVS protein). Different databases use different naming conventions. Without normalization, cross-database matching fails silently.

**Implementation**: `utils/variant_normalizer.py`
- Parse protein change shorthand (V600E → p.Val600Glu) using amino acid code mapping
- Parse common variant formats: cDNA (c.1799T>A), genomic (g.140453136A>T), protein (p.V600E)
- Normalize variant names for cross-database querying (strip whitespace, standardize stop codon notation: * vs Ter vs X)
- Provide HGVS validation for well-formed expressions

**New tool**: `normalize_variant` — converts between variant notation formats

### 10c. Functional Domain Awareness (UniProt/InterPro) 🏗️

**Why**: Oncogenicity code OM1 requires the variant to be "located in critical/functional domain without benign variation." Currently we just check if the gene is a known oncogene/TSG — but the real question is whether amino acid position 600 falls within the BRAF kinase activation segment (it does). This is the difference between "BRAF is important" and "V600 is in the kinase domain."

**Implementation**: `clients/uniprot_client.py`
- Query UniProt REST API (`https://rest.uniprot.org/uniprotkb/search`)
- Fetch protein features: domains, active sites, binding sites, PTM sites, mutagenesis annotations
- Map amino acid position to functional regions
- Check if variant falls in a PFAM/InterPro domain critical for protein function

**New model**: `ProteinDomainInfo` with fields: domain_name, domain_type (kinase, DNA-binding, etc.), start_pos, end_pos, source (UniProt/InterPro), significance

**Integration**: Enhance OM1 auto-detection in oncogenicity scorer with domain-level evidence.

### 10d. PubMed/Entrez Literature Search 📚

**Why**: The roadmap lists NCBI Entrez/PubMed as a data source, but no PubMed client was built. Literature co-occurrence (how many papers mention BRAF AND V600E AND melanoma) is fundamental for variant interpretation. PS4, PP4, and other ACMG criteria reference published studies. CIViC evidence items cite PMIDs — but there may be additional literature not yet curated.

**Implementation**: `clients/pubmed_client.py`
- Use NCBI E-utilities (esearch + efetch) to search PubMed
- Gene + variant co-occurrence search (e.g., `BRAF[gene] AND V600E AND cancer`)
- Return: total publication count, recent papers (title, authors, journal, year, PMID), MeSH terms
- Rate limit: 3 req/sec without key, 10 with NCBI API key (reuse existing env var)

**New tool**: `search_literature` — find publications for a gene/variant/disease combination
**New tool**: `get_publication` — fetch details for a specific PMID

### 10e. In-Silico Prediction Score Aggregation 🤖

**Why**: ACMG criteria PP3 ("multiple computational evidence supports deleterious effect") and BP4 ("computational evidence suggests no impact") require aggregated in-silico prediction scores. Oncogenicity code OP1 ("all computational evidence supports oncogenic effect") and SBP1 ("all computational evidence suggests benign") also need this. These are standard columns in every clinical genomics pipeline.

**Implementation**: `clients/dbnsfp_client.py` or integration via MyVariant.info
- Query MyVariant.info API (`https://myvariant.info/v1/variant/`) which aggregates:
  - **SIFT**: tolerated vs. damaging (score < 0.05 = damaging)
  - **PolyPhen-2**: benign/possibly/probably damaging (score > 0.85 = probably damaging)
  - **REVEL**: meta-predictor (>0.5 likely pathogenic, >0.75 likely disease-causing)
  - **CADD**: Combined Annotation Dependent Depletion (phred > 20 = top 1%)
  - **AlphaMissense**: DeepMind's structure-based pathogenicity (>0.564 likely pathogenic)
  - **GERP++**: Conservation score (>2 = conserved)
  - **phyloP**: Evolutionary conservation
- Return aggregated consensus: how many predictors say deleterious vs. benign

**New model**: `InSilicoPredictions` with per-tool scores and consensus

**New tool**: `predict_variant_effect` — aggregate computational predictions for a variant

### 10f. Enhanced Oncogenicity Auto-Detection 🧮

**Why**: With gnomAD + in-silico + domain data, the auto-detection becomes dramatically more accurate:
- **SBVS1**: gnomAD AF > 5% in any continental population → -8 points (currently not possible)
- **SBS1**: gnomAD AF > 1% → -4 points (currently guessing)
- **OM4**: gnomAD AF < 0.01% → +2 points (currently guessing from ClinVar labels)
- **OP4**: Absent from gnomAD → +1 point (currently impossible)
- **OP1**: All in-silico predictors agree → +1 point (currently not checked)
- **SBP1**: All in-silico predictors benign → -1 point (currently not checked)
- **OM1**: Variant in critical domain (UniProt) → +2 points (currently gene-level only)
- **PP3/BP4** (ACMG): Computational evidence → now has actual scores

### Summary: New Tools Added in Phase 10

| # | Tool | Description |
|---|------|-------------|
| 21 | `lookup_gnomad_frequency` | Population allele frequencies from gnomAD |
| 22 | `normalize_variant` | Convert between variant notation formats (V600E ↔ p.Val600Glu) |
| 23 | `lookup_protein_domains` | Functional domain mapping from UniProt/InterPro |
| 24 | `search_literature` | PubMed literature co-occurrence search |
| 25 | `get_publication` | Fetch publication details by PMID |
| 26 | `predict_variant_effect` | Aggregated in-silico pathogenicity predictions |

### Key References for Phase 10

| Feature | Reference | PMID |
|---------|-----------|------|
| gnomAD v4 | Chen et al. (2024). Nature. | 38862018 |
| HGVS nomenclature | den Dunnen et al. (2016). Hum Mutat. | 26931183 |
| UniProt domains | UniProt Consortium (2023). Nucleic Acids Res. | 36408920 |
| REVEL | Ioannidis et al. (2016). Am J Hum Genet. | 27666373 |
| CADD | Rentzsch et al. (2021). Nucleic Acids Res. | 33237323 |
| AlphaMissense | Cheng et al. (2023). Science. | 37733863 |
| MyVariant.info | Xin et al. (2016). Genome Biol. | 27154141 |
| ClinGen SVI | Pejaver et al. (2022). Am J Hum Genet. | 36413997 |