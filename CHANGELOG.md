# Changelog

All notable changes to CAPVIC are documented in this file.

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
