# 🧬 CAPVIC — Clinical & Academic Precision Variant Interpretation Console

> **A production-grade MCP server for precision oncology variant classification**
> Integrating 8 data sources (CIViC, ClinVar, OncoKB, VICC MetaKB, gnomAD, UniProt, PubMed, MyVariant.info) with AMP/ASCO/CAP, ClinGen/CGC/VICC Oncogenicity SOP, and ACMG/AMP classification frameworks.

[![CI](https://github.com/ArioMoniri/CAPVIC/actions/workflows/ci.yml/badge.svg)](https://github.com/ArioMoniri/CAPVIC/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![MCP](https://img.shields.io/badge/MCP-Compatible-purple.svg)](https://modelcontextprotocol.io)

---

## 📖 Table of Contents

- [Overview](#-overview)
- [Architecture](#-architecture)
- [Quick Start](#-quick-start)
- [Docker](#-docker)
- [Claude Desktop Integration](#-claude-desktop-integration)
- [Environment Variables](#-environment-variables)
- [Tool Reference (26 Tools)](#-tool-reference)
- [Classification Frameworks](#-classification-frameworks)
- [Example Workflows](#-example-workflows)
- [For ML/AI Competition Teams](#-for-mlai-competition-teams)
- [Data Sources & Freshness](#-data-sources--freshness)
- [Development](#-development)
- [Key References](#-key-references)
- [Known Limitations & Future Work](#-known-limitations--future-work)
- [Disclaimer](#%EF%B8%8F-disclaimer)

---

## 🔬 Overview

CAPVIC turns any AI assistant (Claude, GPT, etc.) into a **virtual molecular tumor board** by providing real-time access to the world's leading clinical genomics knowledgebases and implementing gold-standard classification frameworks as computational logic.

### What It Does

| Capability | Description |
|------------|-------------|
| 🔍 **Multi-source evidence search** | Query CIViC, ClinVar, OncoKB, and MetaKB simultaneously |
| 🏷️ **Somatic classification** | AMP/ASCO/CAP 4-tier (Tier I–IV) with evidence levels (A–D) |
| 📊 **Oncogenicity scoring** | ClinGen/CGC/VICC SOP point-based system (18 evidence codes) |
| 🧪 **Germline interpretation** | ACMG/AMP 5-tier framework reference with ClinVar context |
| ⚖️ **Cross-database comparison** | Side-by-side concordance/discordance analysis |
| 📋 **Pathogenicity summary** | Structured pathogenic vs. benign evidence report |
| 💊 **Therapeutic matching** | FDA-approved & investigational therapies from OncoKB/CIViC |
| 🧬 **Population frequency** | gnomAD allele frequencies for BA1/PM2/SBVS1 criteria |
| 🔀 **Variant normalization** | HGVS notation parsing (V600E ↔ p.Val600Glu) |
| 🏗️ **Protein domain mapping** | UniProt/InterPro domain lookup for OM1 evidence |
| 📚 **Literature search** | PubMed co-occurrence for gene/variant/disease |
| 🤖 **In-silico predictions** | SIFT, PolyPhen-2, REVEL, CADD, AlphaMissense aggregation |

### Data Sources

| Source | Type | Access | Update Frequency |
|--------|------|--------|-----------------|
| 🟢 [CIViC](https://civicdb.org) | Expert-curated clinical evidence | Free, no auth | Continuous (community-driven) |
| 🟢 [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Aggregate pathogenicity from 2000+ labs | Free, no auth | Weekly |
| 🟡 [OncoKB](https://www.oncokb.org) | FDA-recognized oncogenicity + therapy levels | Free academic token | Continuously |
| 🟢 [VICC MetaKB](https://search.cancervariants.org) | Harmonized from 6 knowledgebases | Free, no auth | Periodic |
| 🟢 [gnomAD](https://gnomad.broadinstitute.org) | Population allele frequencies (152k genomes) | Free, no auth | Major releases |
| 🟢 [UniProt](https://www.uniprot.org) | Protein domains & functional annotation | Free, no auth | Monthly |
| 🟢 [PubMed](https://pubmed.ncbi.nlm.nih.gov) | Biomedical literature (36M+ articles) | Free, no auth | Daily |
| 🟢 [MyVariant.info](https://myvariant.info) | Aggregated in-silico predictions (dbNSFP) | Free, no auth | Periodic |

---

## 🏗️ Architecture

```
┌────────────────────────────────────────────────────────────────┐
│                        AI Assistant                             │
│              (Claude, GPT, or any MCP client)                  │
└─────────────────────────┬──────────────────────────────────────┘
                          │ MCP (stdio transport)
┌─────────────────────────▼──────────────────────────────────────┐
│                    CAPVIC MCP Server                            │
│                  (FastMCP · Python 3.11+)                       │
│                                                                │
│  ┌──────────────────┐  ┌──────────────────┐  ┌──────────────┐ │
│  │  26 MCP Tools    │  │  Classification   │  │  Formatters  │ │
│  │                  │  │  Engines          │  │              │ │
│  │  • Search        │  │  • AMP/ASCO/CAP   │  │  • Reports   │ │
│  │  • Classify      │  │  • Oncogenicity   │  │  • Tables    │ │
│  │  • Compare       │  │    SOP            │  │  • Evidence  │ │
│  │  • Annotate      │  │  • ACMG/AMP       │  │    summaries │ │
│  │  • Discover      │  │    helper         │  │              │ │
│  │  • Predict       │  │                   │  │              │ │
│  │  • Normalize     │  │                   │  │              │ │
│  └──────┬───────────┘  └──────────────────┘  └──────────────┘ │
│         │                                                      │
│  ┌──────▼───────────────────────────────────────────────────┐  │
│  │              Async API Clients (httpx)                    │  │
│  │   Rate-limited · Retry with backoff · Graceful fallback  │  │
│  └──┬───────┬──────────┬───────────┬─────────┬─────────────┘  │
└─────┼───────┼──────────┼───────────┼─────────┼────────────────┘
      │       │          │           │         │
┌─────▼──┐ ┌─▼────────┐ ┌▼────────┐ ┌▼──────┐ ┌▼───────────────┐
│ CIViC  │ │ ClinVar  │ │ OncoKB │ │MetaKB│ │gnomAD/UniProt/ │
│GraphQL │ │E-utility │ │REST API│ │ REST │ │PubMed/MyVar   │
└────────┘ └──────────┘ └────────┘ └──────┘ └───────────────┘
```

### Project Structure

```
CAPVIC/
├── src/variant_mcp/
│   ├── server.py                  # FastMCP server + 26 tool registrations
│   ├── constants.py               # URLs, enums, evidence codes, disclaimers
│   ├── clients/
│   │   ├── base_client.py         # Async HTTP + rate limiting + retry
│   │   ├── civic_client.py        # CIViC V2 GraphQL client
│   │   ├── clinvar_client.py      # NCBI E-utilities (ClinVar)
│   │   ├── oncokb_client.py       # OncoKB REST API
│   │   ├── metakb_client.py       # VICC MetaKB search
│   │   ├── gnomad_client.py       # gnomAD GraphQL (population frequencies)
│   │   ├── pubmed_client.py       # PubMed E-utilities (literature search)
│   │   ├── uniprot_client.py      # UniProt REST (protein domains)
│   │   └── myvariant_client.py    # MyVariant.info (in-silico predictions)
│   ├── utils/
│   │   └── variant_normalizer.py  # Pure Python HGVS notation parser
│   ├── classification/
│   │   ├── amp_asco_cap.py        # AMP/ASCO/CAP 4-tier somatic classifier
│   │   ├── oncogenicity_sop.py    # ClinGen/CGC/VICC oncogenicity scorer
│   │   └── acmg_amp.py           # ACMG/AMP germline pathogenicity helper
│   ├── models/
│   │   ├── inputs.py             # Pydantic v2 input models (validated)
│   │   ├── evidence.py           # Unified evidence data models
│   │   └── classification.py     # Classification result models
│   ├── formatters/
│   │   ├── reports.py            # Markdown report generators
│   │   └── tables.py            # Framework reference table formatters
│   └── queries/
│       └── civic_graphql.py      # CIViC GraphQL query definitions
├── tests/                        # 114 unit tests (all clients mock-tested)
├── .github/workflows/ci.yml     # CI: lint, typecheck, test (3.11+3.12), build
├── Dockerfile
├── docker-compose.yml
├── pyproject.toml
└── .env.example
```

---

## 🚀 Quick Start

### Prerequisites

- Python 3.11 or later
- (Optional) OncoKB academic API token — [register free](https://www.oncokb.org/apiAccess)
- (Optional) NCBI API key — [get one here](https://www.ncbi.nlm.nih.gov/account/settings/)

### Install

```bash
# Clone the repository
git clone https://github.com/ArioMoniri/CAPVIC.git
cd CAPVIC

# Install in development mode
pip install -e ".[dev]"

# Run the server
python -m variant_mcp.server
```

### Verify Installation

```bash
# Check server loads correctly
python -c "from variant_mcp.server import mcp; print('✅ CAPVIC server OK')"

# Run tests
pytest tests/ -v --tb=short -m "not integration"
```

---

## 🐳 Docker

```bash
# Build and run
docker compose build
docker compose up

# With API tokens
ONCOKB_API_TOKEN=your_token NCBI_API_KEY=your_key docker compose up
```

---

## 🖥️ Claude Desktop Integration

Add to your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "variant_mcp": {
      "command": "python",
      "args": ["-m", "variant_mcp.server"],
      "cwd": "/path/to/CAPVIC",
      "env": {
        "PYTHONPATH": "/path/to/CAPVIC/src",
        "ONCOKB_API_TOKEN": "your_token_here",
        "NCBI_API_KEY": "your_key_here"
      }
    }
  }
}
```

Or using the installed entry point:

```json
{
  "mcpServers": {
    "variant_mcp": {
      "command": "variant-mcp",
      "env": {
        "ONCOKB_API_TOKEN": "your_token_here"
      }
    }
  }
}
```

---

## 🔑 Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `ONCOKB_API_TOKEN` | No | OncoKB API token ([free academic access](https://www.oncokb.org/apiAccess)). Without it, OncoKB tools return a helpful registration message. |
| `NCBI_API_KEY` | No | NCBI API key. Increases ClinVar rate limit from 3 → 10 requests/sec. |

---

## 🛠️ Tool Reference

### Category 1: Multi-Source Power Tools ⚡

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 1 | `variant_search_evidence` | 🔍 Query ALL databases simultaneously for unified evidence | `gene`\*, `variant`, `disease`, `therapy`, `evidence_type`, `sources`, `limit` |
| 2 | `variant_classify` | 🏷️ Apply classification framework(s) — somatic or germline | `gene`\*, `variant`\*, `disease`, `variant_origin`\* |
| 3 | `variant_compare_sources` | ⚖️ Cross-database concordance/discordance analysis | `gene`\*, `variant`\* |
| 20 | `variant_pathogenicity_summary` | 📋 Structured pathogenic vs. benign evidence report | `gene`\*, `variant`\*, `disease` |

### Category 2: CIViC Tools 🧬

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 4 | `civic_search_evidence` | Search CIViC clinical evidence items | `gene`, `variant`, `disease`, `therapy`, `evidence_type` |
| 5 | `civic_search_assertions` | Search curated CIViC assertions | `gene`, `disease`, `therapy`, `significance` |
| 6 | `civic_get_gene` | Get gene details + variant list | `name`\* |
| 7 | `civic_get_variant` | Get variant + molecular profiles | `variant_id`\* |
| 8 | `civic_get_evidence_item` | Get full evidence item detail | `evidence_id`\* |

### Category 3: ClinVar Tools 🏥

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 9 | `clinvar_search` | Search ClinVar by gene, variant, significance, disease | `gene`, `variant`, `clinical_significance` |
| 10 | `clinvar_get_variant` | Full record by variation ID, rsID, or HGVS | `variation_id`, `rsid`, `hgvs` |

### Category 4: OncoKB Tools 💊

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 11 | `oncokb_annotate` | Annotate somatic mutation (oncogenicity + therapy) | `gene`\*, `variant`\*, `tumor_type` |
| 12 | `oncokb_cancer_genes` | List all OncoKB cancer genes (oncogene/TSG status) | — |

### Category 5: Classification Framework Tools 📊

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 13 | `classify_amp_tier` | AMP/ASCO/CAP 4-tier somatic classification | `gene`\*, `variant`\*, `disease` |
| 14 | `score_oncogenicity` | ClinGen/CGC/VICC point-based oncogenicity SOP | `gene`\*, `variant`\*, `evidence_codes` |
| 15 | `explain_acmg_criteria` | ACMG/AMP germline criteria reference | `criteria_code`, `query` |

### Category 6: Discovery Tools 🔎

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 16 | `lookup_gene` | Gene name autocomplete | `query`\* |
| 17 | `lookup_disease` | Disease name autocomplete (DOID, OncoTree) | `query`\* |
| 18 | `lookup_therapy` | Therapy/drug name autocomplete | `query`\* |

### Category 7: Utility Tools 📚

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 19 | `get_classification_frameworks_reference` | Complete reference for all 3 frameworks | — |

### Category 8: Scientific Enhancement Tools 🧬 (Phase 10)

| # | Tool | Description | Key Inputs |
|---|------|-------------|------------|
| 21 | `lookup_gnomad_frequency` | 🧬 Population allele frequencies from gnomAD (BA1/PM2/SBVS1) | `variant_id`\*, `genome_version` |
| 22 | `normalize_variant` | 🔀 Parse & convert variant notation (V600E ↔ p.Val600Glu) | `variant`\*, `gene` |
| 23 | `lookup_protein_domains` | 🏗️ Protein functional domains from UniProt/InterPro (OM1) | `gene`\*, `variant`, `position` |
| 24 | `search_literature` | 📚 PubMed literature co-occurrence search | `gene`\*, `variant`, `disease`, `limit` |
| 25 | `get_publication` | 📄 Fetch full publication details by PMID | `pmid`\* |
| 26 | `predict_variant_effect` | 🤖 SIFT/PolyPhen-2/REVEL/CADD/AlphaMissense predictions (PP3/BP4) | `gene`\*, `variant`\*, `hgvs_id` |

\* = required parameter

---

## 📊 Classification Frameworks

### 1. AMP/ASCO/CAP 4-Tier Somatic Classification

> *Li et al., J Mol Diagn 2017. PMID: [27993330](https://pubmed.ncbi.nlm.nih.gov/27993330/)*

**Use for**: Somatic (cancer) variant clinical significance

| Tier | Name | Evidence Levels | Description |
|------|------|-----------------|-------------|
| **I** | 🔴 Strong Clinical Significance | **A**: FDA-approved, guidelines | **B**: Well-powered studies, consensus |
| **II** | 🟠 Potential Clinical Significance | **C**: FDA-approved (different tumor) | **D**: Preclinical / case reports |
| **III** | 🟡 Unknown Clinical Significance | — | No convincing published evidence |
| **IV** | 🟢 Benign / Likely Benign | — | High population frequency, no cancer association |

### 2. ClinGen/CGC/VICC Oncogenicity SOP

> *Horak et al., Genet Med 2022. PMID: [35101336](https://pubmed.ncbi.nlm.nih.gov/35101336/)*

**Use for**: Somatic variant oncogenicity assessment (point-based)

| Classification | Points | Oncogenic Codes | Benign Codes |
|---------------|--------|-----------------|--------------|
| 🔴 **Oncogenic** | ≥ 10 | OVS1 (+8), OS1-3 (+4 each) | — |
| 🟠 **Likely Oncogenic** | 6–9 | OM1-4 (+2 each) | — |
| 🟡 **VUS** | −5 to +5 | OP1-4 (+1 each) | — |
| 🔵 **Likely Benign** | −6 to −9 | — | SBP1-2 (−1 each) |
| 🟢 **Benign** | ≤ −10 | — | SBVS1 (−8), SBS1-2 (−4 each) |

### 3. ACMG/AMP 5-Tier Germline Pathogenicity

> *Richards et al., Genet Med 2015. PMID: [25741868](https://pubmed.ncbi.nlm.nih.gov/25741868/)*

**Use for**: Germline (inherited) variant pathogenicity

| Classification | Key Criteria |
|---------------|-------------|
| 🔴 **Pathogenic** | PVS1 + ≥1 Strong; or ≥2 Strong; or 1 Strong + ≥3 Moderate/Supporting |
| 🟠 **Likely Pathogenic** | PVS1 + 1 Moderate; or 1 Strong + 1–2 Moderate |
| 🟡 **VUS** | Criteria not met for any other classification |
| 🔵 **Likely Benign** | 1 Strong benign + 1 Supporting benign |
| 🟢 **Benign** | BA1 standalone; or ≥2 Strong benign |

### When to Use Which Framework

| Context | Framework | Notes |
|---------|-----------|-------|
| Somatic variant, clinical actionability question | **AMP/ASCO/CAP** | Focus on therapeutic, diagnostic, prognostic significance |
| Somatic variant, oncogenicity question | **ClinGen/CGC/VICC SOP** | Point-based scoring for oncogenic vs. benign |
| Germline variant, disease causation question | **ACMG/AMP** | For inherited variants — requires expert review |
| Both somatic questions | **AMP + Oncogenicity SOP** | Complementary — AMP = clinical action, SOP = biological mechanism |

---

## 💡 Example Workflows

### 🧪 "Is BRAF V600E pathogenic in melanoma?"

```python
variant_classify(
    gene="BRAF",
    variant="V600E",
    disease="Melanoma",
    variant_origin="somatic"
)
```

Returns: AMP Tier I/A + Oncogenicity Score 14 (Oncogenic) with full evidence trail from CIViC, ClinVar, and OncoKB.

### 💊 "What therapies work for KRAS G12C in NSCLC?"

```python
variant_search_evidence(
    gene="KRAS",
    variant="G12C",
    disease="Non-small Cell Lung Cancer"
)
```

Returns: Sotorasib (Lumakras) and adagrasib (Krazati) from OncoKB Level 1, plus CIViC clinical evidence items.

### 🧬 "Is BRCA1 R1699W pathogenic?"

```python
variant_classify(
    gene="BRCA1",
    variant="R1699W",
    variant_origin="germline"
)
```

Returns: ClinVar aggregate classification with review stars, submitter breakdown, ACMG/AMP criteria context, and conflicting interpretation analysis.

### 📊 "Score TP53 R175H for oncogenicity"

```python
score_oncogenicity(
    gene="TP53",
    variant="R175H"
)
```

Returns: Point breakdown (OVS1: +8, OS3: +4 = 12 points → **Oncogenic**) with evidence code explanations.

### ⚖️ "Compare what all databases say about PIK3CA H1047R"

```python
variant_compare_sources(
    gene="PIK3CA",
    variant="H1047R"
)
```

Returns: Side-by-side concordance across CIViC, ClinVar, OncoKB, and MetaKB showing agreements and discrepancies.

### 🧬 "What is the gnomAD frequency of BRAF V600E?"

```python
lookup_gnomad_frequency(variant_id="7-140453136-A-T")
```

Returns: Global AF, per-population frequencies (7 populations), clinical interpretation (OM4/OP4/BA1/SBS1 applicability). Requires chrom-pos-ref-alt format (see [Limitations](#-known-limitations--future-work)).

### 🤖 "Are computational predictors consistent for TP53 R175H?"

```python
predict_variant_effect(gene="TP53", variant="R175H")
```

Returns: SIFT (0.0, deleterious), PolyPhen-2 (1.0, probably_damaging), REVEL (0.98), CADD (35), AlphaMissense (0.99) → **Consensus: DAMAGING** → supports PP3 + OP1.

### 🏗️ "Is BRAF V600E in a functional domain?"

```python
lookup_protein_domains(gene="BRAF", variant="V600E")
```

Returns: Position 600 falls in the Protein Kinase domain (457–717, Pfam) → supports **OM1** (+2 points in oncogenicity SOP).

### 📚 "Find recent papers on EGFR T790M resistance"

```python
search_literature(gene="EGFR", variant="T790M", disease="lung cancer", limit=5)
```

Returns: Total publications, recent papers with titles/authors/PMIDs for osimertinib resistance.

---

## 🏆 For ML/AI Competition Teams

CAPVIC is designed to help teams building ML models that predict variant pathogenicity (e.g., Kaggle competitions, CAGI challenges). Here's how to use it:

### Validate Predictions

```
Your model predicts EGFR L858R = Pathogenic (score 0.97)

→ Use variant_pathogenicity_summary(gene="EGFR", variant="L858R")
→ CAPVIC returns: ClinVar 30/30 Pathogenic (★★★), OncoKB Oncogenic Level 1,
  CIViC 25+ evidence items — STRONG agreement with your prediction ✅
```

### Find Edge Cases

```
Your model predicts VHL R167Q = VUS (score 0.52)

→ Use variant_classify(gene="VHL", variant="R167Q", variant_origin="germline")
→ CAPVIC returns: ClinVar has conflicting interpretations (3 Pathogenic, 2 VUS),
  review status ★ — genuinely uncertain, your model captures the ambiguity ✅
```

### Understand Failures

```
Your model predicts IDH1 R132H = Benign (score 0.15) — WRONG!

→ Use variant_compare_sources(gene="IDH1", variant="R132H")
→ CAPVIC shows: Every source says Oncogenic/Pathogenic with high confidence.
  Hotspot mutation, gain-of-function → investigate feature engineering 🔍
```

---

## 📡 Data Sources & Freshness

| Source | Update Frequency | Coverage | Notes |
|--------|-----------------|----------|-------|
| **CIViC** | Continuous | ~3,800 evidence items, ~400 genes | Community-curated, expert-reviewed |
| **ClinVar** | Weekly (Sunday) | 2.5M+ variant submissions | Aggregate from 2,000+ laboratories |
| **OncoKB** | Continuously | ~700 genes, ~5,400 alterations | FDA-recognized, MSK-curated |
| **VICC MetaKB** | Periodic | Harmonized from 6 KBs | CIViC + OncoKB + CGI + JAX-CKB + MMatch + PMKB |
| **gnomAD** | Major releases | 76,156+ genomes (v4.1) | Population AF for BA1/PM2/SBVS1/OM4 criteria |
| **UniProt** | Monthly | 570K+ reviewed entries | Protein domains, functional annotation |
| **PubMed** | Daily | 36M+ biomedical articles | Literature co-occurrence for evidence support |
| **MyVariant.info** | Periodic | Aggregated dbNSFP | SIFT, PolyPhen, REVEL, CADD, AlphaMissense |

> 📌 All API responses include retrieval timestamps. ClinVar data may be up to 7 days old; CIViC and OncoKB are near real-time.

---

## 🧑‍💻 Development

### Setup

```bash
pip install -e ".[dev]"
```

### Commands

```bash
# Run tests (114 unit tests)
pytest tests/ -v --tb=short -m "not integration"

# Lint
ruff check src/ tests/
ruff format src/ tests/

# Type check
mypy src/variant_mcp/ --ignore-missing-imports

# Build package
python -m build
```

### CI/CD

GitHub Actions runs on every push:
1. **Lint** — ruff check + format verification
2. **Typecheck** — mypy strict(ish) checking
3. **Test** — pytest on Python 3.11 + 3.12 matrix
4. **Build** — package build + wheel verification

---

## 📚 Key References

| # | Reference | PMID |
|---|-----------|------|
| 1 | Li MM, et al. (2017). "Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer." *J Mol Diagn*, 19(1):4-23. | [27993330](https://pubmed.ncbi.nlm.nih.gov/27993330/) |
| 2 | Horak P, et al. (2022). "Standards for the classification of pathogenicity of somatic variants in cancer (oncogenicity)." *Genet Med*, 24(4):986-998. | [35101336](https://pubmed.ncbi.nlm.nih.gov/35101336/) |
| 3 | Richards S, et al. (2015). "Standards and guidelines for the interpretation of sequence variants." *Genet Med*, 17(5):405-424. | [25741868](https://pubmed.ncbi.nlm.nih.gov/25741868/) |
| 4 | Griffith M, et al. (2017). "CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer." *Nat Genet*, 49(2):170-174. | [28138153](https://pubmed.ncbi.nlm.nih.gov/28138153/) |
| 5 | Wagner AH, et al. (2020). "A harmonized meta-knowledgebase of clinical interpretations of somatic genomic variants in cancer." *Nat Genet*, 52(4):448-457. | [32246132](https://pubmed.ncbi.nlm.nih.gov/32246132/) |
| 6 | Chakravarty D, et al. (2017). "OncoKB: A Precision Oncology Knowledge Base." *JCO Precis Oncol*, 2017:PO.17.00011. | [28890946](https://pubmed.ncbi.nlm.nih.gov/28890946/) |
| 7 | Chen S, et al. (2024). "A genomic mutational constraint map using variation in 76,156 human genomes." *Nature*, 625:92-100. | [38862018](https://pubmed.ncbi.nlm.nih.gov/38862018/) |
| 8 | Ioannidis NM, et al. (2016). "REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants." *Am J Hum Genet*, 99(4):877-885. | [27666373](https://pubmed.ncbi.nlm.nih.gov/27666373/) |
| 9 | Cheng J, et al. (2023). "Accurate proteome-wide missense variant effect prediction with AlphaMissense." *Science*, 381(6664):eadg7492. | [37733863](https://pubmed.ncbi.nlm.nih.gov/37733863/) |
| 10 | Pejaver V, et al. (2022). "Calibration of computational tools for missense variant pathogenicity classification." *Am J Hum Genet*, 109(12):2163-2177. | [36413997](https://pubmed.ncbi.nlm.nih.gov/36413997/) |

---

## 🚧 Known Limitations & Future Work

The following integrations were evaluated but **intentionally excluded** or scoped down to maintain production-readiness. Each would require substantial new code, external service dependencies, or complex parsing that cannot be reliably unit-tested:

| Limitation | Reason | Impact | Workaround |
|-----------|--------|--------|------------|
| **gnomAD gene+variant search** | gnomAD's search API returns HTML URLs that require fragile string parsing to extract variant IDs. No stable API for protein-change → genomic-coordinate mapping exists. | Users must provide `variant_id` in chrom-pos-ref-alt format (e.g., `7-140453136-A-T`) rather than gene+protein change. | Use ClinVar's HGVS expressions or external tools like [Ensembl VEP](https://www.ensembl.org/vep) to convert protein changes to genomic coordinates. |
| **Automated HGVS → genomic coordinate mapping** | Requires a reference genome + transcript database (e.g., UTA, SeqRepo) — heavy infrastructure not suitable for an MCP server. Libraries like `hgvs` (biocommons) need PostgreSQL + 20GB+ of sequence data. | Cannot auto-convert `p.Val600Glu` to `chr7:g.140453136A>T` for gnomAD/MyVariant.info lookups. | Provide genomic HGVS IDs directly, or use the variant normalizer for protein-level notation and then look up genomic coordinates externally. |
| **ClinGen Allele Registry integration** | Would provide canonical allele IDs for cross-database linking. Requires complex registration workflows and the API has no stable versioning. | No canonical allele ID cross-linking between databases. | Use gene+variant name search (works for CIViC, ClinVar, OncoKB) or provide database-specific IDs. |
| **Functional assay databases (MAVE, DMS)** | Multiplexed Assay of Variant Effect data (e.g., from MaveDB) would strengthen OS2 evidence code. API is experimental and data format varies per assay. | No direct functional assay score integration for OS2 evidence. | Cite functional studies found via `search_literature` tool. |
| **SpliceAI / splice prediction** | Would strengthen PVS1 and splice-site variant assessment. Requires either a 10GB+ model or a paid API, and splice prediction is computationally intensive. | Splice-site variants are detected but not scored for splicing impact. | Use external splice predictors (SpliceAI, MaxEntScan) and feed results into oncogenicity assessment manually. |
| **Automated ACMG/AMP germline classification scoring** | Full ACMG/AMP requires combining 28 evidence criteria with complex combination rules (PVS1 decision trees, segregation data, de novo confirmation). Building a validated, rule-complete implementation is a multi-month effort. | ACMG criteria are explained and referenced, but not auto-scored. Germline classification relies on ClinVar's aggregate expert consensus. | Use `explain_acmg_criteria` for criteria reference, and `clinvar_search` for expert-reviewed germline classifications. |
| **PharmGKB / CPIC pharmacogenomics** | Would add drug-gene interaction annotations. Different data model from somatic classification; would require a separate classification framework. | No pharmacogenomic annotations for germline drug metabolism variants. | Use OncoKB for somatic therapeutic levels, which covers FDA-approved biomarkers. |

### 🔬 Bioinformatician's Assessment

From a clinical genomics perspective, the current tool covers the **core evidence gathering and classification workflow** used in molecular tumor boards:

1. **Evidence aggregation** ✅ — CIViC, ClinVar, OncoKB, MetaKB provide the primary sources used in clinical reporting
2. **Somatic classification** ✅ — AMP/ASCO/CAP + oncogenicity SOP cover the standard-of-care frameworks
3. **Population frequency** ✅ — gnomAD direct lookup handles the most common BA1/PM2 queries
4. **In-silico predictions** ✅ — 7 predictors via MyVariant.info covers PP3/BP4 requirements
5. **Protein domain context** ✅ — UniProt provides OM1 evidence
6. **Literature context** ✅ — PubMed search supports evidence review

The limitations above represent **advanced features** that typically require institutional infrastructure (reference genomes, local databases, licensed tools) and are beyond the scope of a lightweight MCP server.

---

## ⚠️ Disclaimer

> **RESEARCH DISCLAIMER**: This tool aggregates data from CIViC, ClinVar, OncoKB, and VICC MetaKB for **research and educational purposes only**. Classification suggestions follow published AMP/ASCO/CAP, ClinGen/CGC/VICC, and ACMG/AMP frameworks but are **computationally derived approximations** — NOT validated clinical interpretations. Clinical variant classification MUST be performed by qualified molecular pathologists using validated laboratory procedures. Do not use for clinical decision-making without professional review.

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

---

<p align="center">
  Built with 🧬 for the precision oncology research community
</p>
