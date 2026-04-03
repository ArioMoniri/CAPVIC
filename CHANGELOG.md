# Changelog

All notable changes to CAPVIC are documented in this file.

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
