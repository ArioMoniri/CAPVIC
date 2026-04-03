# Changelog

All notable changes to CAPVIC are documented in this file.

## [1.0.6] — 2026-04-03

### Added
- **Cancer Hotspots API client** — queries cancerhotspots.org (Chang et al. 2016, 2018) for recurrent somatic mutation hotspots. Uses POST `/api/hotspots/single/byGene` with JSON array body. Returns residue-level recurrence data with sample counts, variant amino acid distributions, cancer type breakdowns, and statistical q-values
- **LitVar2 API client** — queries NCBI LitVar2 (Allot et al. 2023) for variant-centric literature mining. Chains autocomplete search → variant detail → v1 entity enrichment to provide publication counts, disease co-mentions, clinical significance, ClinGen IDs, and all observed HGVS nomenclatures
- **`lookup_cancer_hotspots` tool** — query Cancer Hotspots database by gene or specific residue position. Returns hotspot table with sample counts, q-values, and cancer type distributions
- **`assess_driver_mutation` tool** — cross-references 7 evidence signals (Cancer Hotspots recurrence, OncoKB oncogenic status, CIViC evidence volume, oncogenicity SOP score, in-silico predictions, gnomAD population frequency, known gene role) to produce a composite driver/passenger classification with confidence score (Driver ≥0.7, Likely Driver 0.4–0.69, VUS 0.2–0.39, Passenger <0.2)
- **`search_variant_literature_litvar` tool** — NLP-curated variant literature search via LitVar2. Returns publication count, top disease co-mentions with frequencies, HGVS nomenclature variants found in literature, rsID mapping, and ClinGen IDs
- **`CancerHotspot` model** — Pydantic model for hotspot residue data with sample counts, amino acid distributions, cancer type breakdowns, q-values
- **`DriverMutationAssessment` model** — composite driver classification result with 7-signal scoring, therapeutic implications, and confidence levels
- **Cancer Hotspots integration in evidence bundle** — `_gather_evidence()` now fetches hotspot data in parallel with other sources
- **Enhanced oncogenicity scorer** — `_is_hotspot()` now uses real Cancer Hotspots API data (sample count ≥10 or q-value <0.05) as primary signal, with CIViC/OncoKB as fallbacks
- **53 new unit tests** — comprehensive mock-based tests for Cancer Hotspots client (POST API, parsing, edge cases), LitVar2 client (autocomplete, entity, full search), driver assessment models, hotspot-aware oncogenicity scoring, and GRCh37/38 coordinate awareness

### Changed
- **10 data sources** (up from 8): CIViC, ClinVar, OncoKB, MetaKB, gnomAD, UniProt, PubMed, MyVariant.info, **Cancer Hotspots**, **LitVar2**
- **29 MCP tools** (up from 26): added `lookup_cancer_hotspots`, `assess_driver_mutation`, `search_variant_literature_litvar`
- **Base HTTP client** — widened `json_body` type from `dict` to `Any` to support Cancer Hotspots POST-with-array API pattern

### Fixed
- **Cancer Hotspots API method** — corrected from GET to POST (API requires `POST /hotspots/single/byGene` with JSON array body `["GENE"]`)
- **Cancer Hotspots q-value parsing** — API returns q-values as strings (e.g., `"3.65e-82"`) or `null`; parser now correctly converts to `float | None`
- **LitVar2 `variant/get` trailing slash** — endpoint returns 404 with trailing `/`; removed from URL path
- **LitVar2 `##` URL encoding** — LitVar IDs contain `##` which httpx treats as URL fragment separators, silently stripping them; now pre-encoded as `%23%23`

## [1.0.5] — 2026-04-03

### Fixed
- **CIViC `civic_get_gene` BROKEN**: CIViC deprecated `genes(name:)` query — migrated to `gene(entrezSymbol:)` and `fullName` field (was `officialName`)
- **CIViC `lookup_gene` BROKEN**: `geneTypeahead` field removed from schema — migrated to `featureTypeahead(queryTerm, featureType: GENE)`
- **CIViC `civic_search_assertions` BROKEN**: `significance` parameter changed from `String` to `AssertionSignificance` enum — client now auto-uppercases values
- **ClinVar efetch**: Migrated to VCV format (`rettype=vcv`) after NCBI deprecated `rettype=variation`
- **Test fix**: Updated typeahead test mock from `geneTypeahead` to `featureTypeahead`

## [1.0.4] — 2026-04-03

### Added
- **`output_format` parameter**: All 26 tools now accept `output_format` — `"markdown"` (default, rich formatted), `"json"` (raw structured data from Pydantic models), or `"text"` (plain text, no formatting). Enables programmatic access and pipeline integration
- **Disease/therapy NLP aliases**: `variant_search_evidence`, `civic_search_evidence`, `civic_search_assertions`, `classify_amp_tier`, and `variant_pathogenicity_summary` now auto-expand common abbreviations — CRC → Colorectal Cancer, NSCLC → Lung Non-small Cell Carcinoma, keytruda → Pembrolizumab, gleevec → Imatinib, etc. (35+ disease aliases, 40+ therapy aliases including brand names)

## [1.0.3] — 2026-04-03

### Changed
- **Tool descriptions**: All 26 MCP tool docstrings rewritten for natural language discovery — each now includes a rich description, "Use for" example prompts, and structured output notes. Tools are discoverable via conversational queries like "Find evidences for colorectal cancer therapies involving KRAS mutations"
- **OpenCode config**: Added `ONCOKB_API_TOKEN` and `NCBI_API_KEY` environment variables with `{env:VAR_NAME}` substitution syntax

### Added
- **Natural Language Prompts** section in README with 10 example prompts showing how AI clients map questions to tools
- **Visualization & Artifacts** guide documenting structured markdown outputs (tables, score breakdowns, domain maps) suitable for chart/plot generation by AI clients

## [1.0.2] — 2026-04-03

### Fixed
- **CRITICAL — Oncogenicity thresholds**: Corrected Benign from ≤-10 to ≤-7 and Likely Benign from ≤-6 to ≤-1 per Horak et al. 2022 Table 3 (PMID: 35101336). Previous thresholds created a dead zone where scores -7 to -9 were "Likely Benign" instead of "Benign"
- **REVEL thresholds**: Corrected to exact Pejaver et al. 2022 calibration — PP3_strong ≥0.932 (was ≥0.773), PP3_moderate ≥0.773 (was ≥0.644), PP3_supporting ≥0.644 (was ≥0.5), BP4_strong ≤0.016 (was ≤0.183), BP4_moderate ≤0.183 (was ≤0.290), BP4_supporting ≤0.290 (was ≤0.4)
- **AMP tier parser**: Fixed prefix collision where "TIERI" matched before "TIERII" — now checks longest tier strings first (IV, III, II, I)
- **ClinVar Tier IV**: Conflicting "pathogenic/likely benign" classifications no longer incorrectly assigned Tier IV
- **Rate limiter**: Replaced deprecated `asyncio.get_event_loop()` with `asyncio.get_running_loop()`
- **ClinVar efetch**: NCBI deprecated `rettype=variation` — migrated to `rettype=vcv` with `is_variationid=true`. Rewrote XML parser for VCV format (ClassifiedRecord, structured HGVS sub-elements, Classifications block)
- **Reports assert**: Replaced bare `assert` statements with proper None checks (safe under Python -O flag)
- **Null variant types**: Added `stop_gained` (Sequence Ontology term) and `initiator_codon_variant` to null variant detection

### Added
- **OpenCode integration**: Full configuration guide for [OpenCode](https://opencode.ai/) MCP support with local stdio transport, Docker, and CLI verification commands

## [1.0.1] — 2026-04-03

### Fixed
- **ClinVar esummary parser**: Migrated to new ClinVar API response structure — `clinical_significance` field is now deprecated, data moved to `germline_classification`, `oncogenicity_classification`, and `clinical_impact_classification` fields. Parser now reads all three and combines them.
- **BaseClient.post()**: Was sending GET instead of POST (critical HTTP method bug)
- **CIViC gene extraction**: Gene name was parsed from variant name instead of molecularProfile name (returned "V600E" as gene instead of "BRAF")
- **CIViC GraphQL errors**: Error messages now extract individual message strings instead of dumping raw list
- **ClinVar rsID search**: Fixed field tag from `[dbsnp_id]` to `[dbsnp id]` (NCBI uses space, not underscore)
- **ClinVar efetch**: Added explicit `retmode=xml` parameter
- **OncoKB HGVSg**: Default genome changed from GRCh37 to GRCh38; added try/except error handling
- **MetaKB exception handling**: Narrowed from broad `Exception` catch to specific types
- **Oncogenicity SOP**: OM4/OP4/SBVS1/SBS1 evidence codes now use actual gnomAD allele frequency data instead of incorrectly inferring from ClinVar classification (Horak et al. 2022)
- **Oncogenicity SOP**: Added SBVS1 auto-detection (MAF >5% threshold) — was previously missing entirely
- **Oncogenicity SOP**: Removed arbitrary `len(applied) < 4` limit on OP2 code application
- **gnomAD AF computation**: Added explicit `an > 0` division-by-zero guards
- **Constants**: Removed incorrect "2025 update" reference from AMP/ASCO/CAP comment
- **ClinVar review stars**: Added missing mappings for "criteria provided, conflicting classifications" (★), "no classification provided" (☆), and "no classifications from unflagged records" (☆)
- **MetaKB endpoint**: Fixed API endpoint from non-existent `api/v2/search` to correct `api/v1/associations` (GA4GH spec)
- **MetaKB parser**: Updated response parser for v1 format (`hits.hits[].association`) with correct field extraction for phenotypes, environmentalContexts, evidence source, and evidence_label

### Added
- `gnomad_frequency` field and `has_gnomad_data` property on EvidenceBundle model
- gnomAD source listed in `sources_queried` when frequency data is present
- Classification report now suggests fetching gnomAD data when missing for better SBVS1/SBS1/OM4/OP4 scoring
- Pathogenicity summary notes when gnomAD data is unavailable
- Genome build documentation in tool docstrings (gnomAD, MyVariant.info, OncoKB)
- Updated Genome Build Reference table in README with ClinVar and tips

### Changed
- ROADMAP.md rewritten from build prompt to structured completed/planned roadmap
- README updated with genome build guidance and classification workflow tips

## [1.0.0] — 2026-04-03

### Added
- 26 MCP tools for variant classification, evidence search, and annotation
- 8 data source integrations: CIViC, ClinVar, OncoKB, VICC MetaKB, gnomAD, UniProt, PubMed, MyVariant.info
- AMP/ASCO/CAP 4-tier somatic classification engine
- ClinGen/CGC/VICC oncogenicity SOP scorer (18 evidence codes, point-based)
- ACMG/AMP germline pathogenicity reference helper
- Variant normalizer (HGVS parser: V600E <-> p.Val600Glu)
- UniProt protein domain lookup with OM1 auto-detection
- In-silico prediction aggregation (SIFT, PolyPhen-2, REVEL, CADD, AlphaMissense, GERP, phyloP)
- ClinGen SVI-calibrated REVEL thresholds for PP3/BP4 evidence strength (Pejaver et al. 2022)
- PubMed literature co-occurrence search
- gnomAD population frequency lookup (GRCh37 + GRCh38)
- Markdown report formatters (evidence reports, pathogenicity summaries, source comparisons)
- Evidence bundle with integrated domain + prediction data fed into classifiers
- Docker support with stdio transport
- Claude Desktop integration configs (Docker + native)
- GitHub Actions CI/CD (lint, typecheck, test matrix 3.11+3.12, build)
- 131 unit tests with full mock coverage

### Fixed
- CIViC GraphQL API: migrated from deprecated `geneName`/`variantName` to `molecularProfileName` filter
- gnomAD GraphQL API: removed deprecated `af` field, compute AF from `ac/an`
- UniProt API: fixed `source` field parsing (string vs dict), added `reviewed:true` filter for canonical entries
- MyVariant.info: fixed search query to use structured `dbnsfp.aa.pos`/`ref`/`alt` fields
- Dockerfile: fixed build order (copy source before pip install)

### References
- Li MM et al. (2017) AMP/ASCO/CAP. J Mol Diagn. PMID: 27993330
- Horak P et al. (2022) ClinGen/CGC/VICC oncogenicity SOP. Genet Med. PMID: 35101336
- Richards S et al. (2015) ACMG/AMP. Genet Med. PMID: 25741868
- Pejaver V et al. (2022) ClinGen SVI REVEL calibration. Am J Hum Genet. PMID: 36413997
