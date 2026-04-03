# CAPVIC Roadmap

## Completed (v1.0.5 — 2026-04-03)

### CIViC API Migration Fixes (v1.0.4 → v1.0.5)
- [x] Fixed `civic_get_gene` — migrated from deprecated `genes(name:)` to `gene(entrezSymbol:)`
- [x] Fixed `lookup_gene` — migrated from removed `geneTypeahead` to `featureTypeahead`
- [x] Fixed `civic_search_assertions` — `significance` now uses `AssertionSignificance` enum type
- [x] Fixed ClinVar efetch VCV format migration (`rettype=vcv`)
- [x] Updated test mocks for new CIViC GraphQL schema

## Completed (v1.0.4 — 2026-04-03)

### Output Format & NLP Enhancement (v1.0.3 → v1.0.4)
- [x] Added `output_format` parameter to all 26 tools (markdown/json/text)
- [x] JSON format returns raw Pydantic model data for programmatic access
- [x] Text format strips markdown formatting for plain text consumers
- [x] Added disease alias normalization (35+ abbreviations: CRC, NSCLC, GBM, TNBC, etc.)
- [x] Added therapy alias normalization (40+ brand-to-generic mappings: keytruda → Pembrolizumab, etc.)
- [x] Enhanced variant_search_evidence, civic_search_evidence, civic_search_assertions with NLP normalization

## Completed (v1.0.3 — 2026-04-03)

### Natural Language & Discovery UX (v1.0.2 → v1.0.3)
- [x] Enhanced all 26 tool descriptions for natural language discovery (nexonco-style)
- [x] Added "Use for" prompt examples to every tool docstring
- [x] Added structured markdown output notes for visualization/artifact generation
- [x] Added API key environment variables to OpenCode config (ONCOKB_API_TOKEN, NCBI_API_KEY)
- [x] Added Natural Language Prompts section to README with 10 example prompts
- [x] Added Visualization & Artifacts guide showing structured output capabilities

## Completed (v1.0.2 — 2026-04-03)

### Bug Fixes & Scientific Accuracy (v1.0.1 → v1.0.2)
- [x] Corrected oncogenicity SOP thresholds to match Horak et al. 2022 Table 3
- [x] Corrected REVEL PP3/BP4 thresholds to match Pejaver et al. 2022 exact calibration
- [x] Fixed AMP tier parser prefix collision (TIERI/TIERII)
- [x] Fixed ClinVar API migration (3 new classification fields)
- [x] Fixed MetaKB endpoint and response parser for v1 API
- [x] Fixed CIViC gene extraction from molecularProfile
- [x] Added gnomAD-based SBVS1/SBS1/OM4/OP4 evidence codes
- [x] Added missing null variant types (stop_gained, initiator_codon_variant)

## Completed (v1.0.0 — 2026-04-03)

### Core Platform
- [x] FastMCP server with stdio transport (26 MCP tools)
- [x] Docker support with Claude Desktop integration
- [x] GitHub Actions CI/CD (lint, typecheck, test matrix Python 3.11 + 3.12)
- [x] 131 unit tests with full mock coverage

### Data Source Integrations (8 sources)
- [x] **CIViC** — GraphQL v2 API (molecularProfileName filter, evidence + assertions)
- [x] **ClinVar** — NCBI E-utilities (esearch, esummary, efetch with XML parsing)
- [x] **OncoKB** — REST API with graceful degradation when no token
- [x] **VICC MetaKB** — Harmonized search across 6 knowledgebases
- [x] **gnomAD** — GraphQL v4 (GRCh38) / v2 (GRCh37), AF computed from ac/an
- [x] **UniProt** — Protein domain lookup with reviewed:true canonical filtering
- [x] **PubMed** — Literature co-occurrence search via E-utilities
- [x] **MyVariant.info** — In-silico prediction aggregation (SIFT, PolyPhen-2, REVEL, CADD, AlphaMissense, GERP, phyloP)

### Classification Frameworks
- [x] **AMP/ASCO/CAP 4-Tier** somatic classification (Li et al. 2017)
- [x] **ClinGen/CGC/VICC Oncogenicity SOP** point-based scorer (Horak et al. 2022, 18 evidence codes)
- [x] **ACMG/AMP 5-Tier** germline pathogenicity reference helper (Richards et al. 2015)
- [x] ClinGen SVI-calibrated REVEL thresholds for PP3/BP4 (Pejaver et al. 2022)

### Utilities
- [x] Variant normalizer (HGVS parser: V600E <-> p.Val600Glu)
- [x] UniProt protein domain lookup with OM1 auto-detection
- [x] In-silico prediction aggregation with PP3/BP4 calibrated strength
- [x] Markdown report formatters (evidence reports, pathogenicity summaries, source comparisons)
- [x] Evidence bundle pipeline — domain + prediction data fed into classifiers

---

## Planned (v1.1.0)

### New Data Sources
- [ ] **COSMIC** — Catalogue of Somatic Mutations in Cancer (hotspot frequencies, mutational signatures)
- [ ] **PharmGKB** — Pharmacogenomics associations for drug-variant interactions
- [ ] **dbSNP** — rsID-based variant lookups and cross-referencing
- [ ] **InterVar** — Automated ACMG/AMP interpretation

### Classification Enhancements
- [ ] Full ACMG/AMP 5-tier automated classifier (beyond reference helper)
- [ ] AMP/ASCO/CAP 2025 update integration (Chakravarty et al.)
- [ ] SpliceAI integration for splice-affecting variant detection (PVS1 support)
- [ ] Gene-specific ACMG criteria (ClinGen Variant Curation Expert Panels)

### Platform Features
- [ ] HTTP/SSE transport mode (for web deployment, multi-client)
- [ ] Batch variant processing (VCF file input)
- [ ] Result caching layer (reduce redundant API calls)
- [ ] Structured JSON output mode (for ML pipeline integration)
- [ ] WebSocket real-time updates for long-running batch jobs

---

## Planned (v1.2.0)

### ML/AI Integration
- [ ] Feature vector export (fixed-width numeric arrays for model training)
- [ ] VCF-to-feature-matrix batch pipeline
- [ ] Pre-computed training dataset generation (ClinVar truth set)
- [ ] Model evaluation harness (concordance with expert classifications)

### Clinical Workflow
- [ ] Multi-variant report generation (panel/exome results)
- [ ] Variant-of-the-day training tool for molecular pathology education
- [ ] Audit trail / provenance tracking for classification decisions
- [ ] FHIR-compatible output for EHR integration

### Infrastructure
- [ ] Helm chart for Kubernetes deployment
- [ ] OpenTelemetry observability (latency, error rates, API availability)
- [ ] Rate limit dashboard and quota management
- [ ] Automated API health monitoring with fallback routing

---

## Contributing

See [CHANGELOG.md](CHANGELOG.md) for release history. Contributions welcome — open an issue to discuss new data sources or classification features before submitting a PR.
