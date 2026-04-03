# Example Prompts with Real API Answers

Real responses from CAPVIC v1.0.6 tools against live APIs (tested 2026-04-03).

---

## LitVar2 Literature Mining

### Prompt: "How many papers mention BRAF V600E? What diseases is it associated with?"

**Tool**: `search_variant_literature_litvar(gene="BRAF", variant="V600E")`

**Real answer**:
```
found: True
rsid: rs113488022
pmids_count: 31630
clinical_significance: [other, drug-response, pathogenic, likely-pathogenic]
clingen_ids: [CA16602736, CA281998, CA123643]
chromosome_position: [7:140753336]  (GRCh38)
snp_class: [snv]

Top disease co-mentions (from 20,572 NLP-mined publications):
  Melanoma: 26,551
  Neoplasms: 18,463
  Colorectal Neoplasms: 6,505
  Thyroid Neoplasms: 3,941
  Thyroid cancer, papillary: 3,710

HGVS nomenclatures found in literature:
  p.V600E: 20,322 mentions
  c.1799T>A: 1,025 mentions
  p.V600E,K: 280 mentions

First published: 2004
Publication years span: 2004-2022
```

### Prompt: "Find literature on TP53 R175H mutations"

**Tool**: `search_variant_literature_litvar(gene="TP53", variant="R175H")`

**Real answer**:
```
found: True
rsid: rs28934578
pmids_count: 3,061
clinical_significance: [likely-pathogenic, pathogenic, not-provided]

Top disease co-mentions:
  Neoplasms: 1,150
  Breast Neoplasms: 223
  Ovarian Neoplasms: 154
  Adenomatous Polyposis Coli: 114
  Lung Neoplasms: 107

HGVS nomenclatures:
  p.R175H: 1,362 mentions
  c.524G>A: 97 mentions

First published: 1989
```

### Prompt: "What's in the literature about KRAS G12D?"

**Tool**: `search_variant_literature_litvar(gene="KRAS", variant="G12D")`

**Real answer**:
```
found: True
rsid: rs121913529
pmids_count: 15,089
clinical_significance: [pathogenic, not-provided, likely-pathogenic]

Top disease co-mentions:
  Neoplasms: 5,211
  Pancreatic Neoplasms: 1,524
  Lung Neoplasms: 1,397
  Colorectal Neoplasms: 1,204
  Carcinoma, Non-Small-Cell Lung: 674

HGVS nomenclatures:
  p.G12D: 5,892 mentions
  p.G12V: 2,685 mentions (note: same rsID covers multiple G12 substitutions)
  c.12G>A: 614 mentions

First published: 1987
```

### Prompt: "Search for a variant that doesn't exist"

**Tool**: `search_variant_literature_litvar(gene="FAKEGENE", variant="X999Z")`

**Real answer**:
```
found: False
gene: FAKEGENE
variant: X999Z
pmids_count: 0
```

---

## Cancer Hotspots

### Prompt: "What are the mutation hotspots in BRAF?"

**Tool**: `lookup_cancer_hotspots(gene="BRAF")`

**Real answer** (12 hotspot residues found):
```
Top hotspots:
  V600: 897 samples, q=0.00e+00
    Variant amino acids: E (833), M (29), K (24), R (4), V (4), G (3)
    Top cancer types: skin (357), thyroid (316), bowel (113), lung (33)

  G469: 61 samples, q=3.65e-82
    Variant amino acids: A (34), V (11), R (8), E (4), S (2)
    Top cancer types: lung (23), skin (13), bowel (6)

  K601: 38 samples, q=1.95e-70
    Variant amino acids: E (32), N (5), T (1)
    Top cancer types: skin (8), lung (6), prostate (6)

  D594: 29 samples, q=6.72e-56
  N581: 12 samples, q=2.87e-26
  G466: 11 samples, q=7.67e-18
  ... (12 total)
```

### Prompt: "Is BRAF V600 a known cancer hotspot?"

**Tool**: `lookup_cancer_hotspots(gene="BRAF", residue=600)`

**Real answer**:
```
gene: BRAF
residue: V600
total_sample_count: 897
classification: single residue
q_value: 0.0  (highly significant)
q_value_cancertype: 0.0

Variant amino acid distribution:
  E: 833 (92.9%)
  M: 29 (3.2%)
  K: 24 (2.7%)
  R: 4, V: 4, G: 3

Cancer type distribution:
  skin: 357
  thyroid: 316
  bowel: 113
  lung: 33
  brain: 20
  unk: 21
  ... (18 cancer types total)
```

### Prompt: "What are the KRAS hotspots?"

**Tool**: `lookup_cancer_hotspots(gene="KRAS")`

**Real answer** (12 hotspot residues found):
```
  G12: 2,175 samples, q=0.00e+00
    Amino acids: D (757), V (657), C (384), R (166), A (134), S (68)
    Top types: pancreas (745), lung (570), bowel (500)

  G13: 264 samples, q=0.00e+00
    Amino acids: D (214), C (32), V (6), G (6)
    Top types: bowel (130), lung (57)

  Q61: 190 samples, q=0.00e+00
    Amino acids: H (108), R (37), L (21), K (20)
    Top types: pancreas (64), bowel (31), lung (28)

  A146: 61 samples
  K117: 13 samples
  ... (12 total)
```

### Prompt: "Show TP53 hotspots"

**Tool**: `lookup_cancer_hotspots(gene="TP53")`

**Real answer** (120 hotspot residues found — TP53 is extensively mutated):
```
Top 5:
  R273: 609 samples — brain (113), bowel (105), lung (78)
  R248: 560 samples — bowel (98), stomach (70), brain (54)
  R175: 416 samples — bowel (108), stomach (63), pancreas (39)
  Y220: 256 samples — bowel (76), lung (26), stomach (21)
  R249: 230 samples — liver (138), stomach (20), brain (8)

Note: TP53 has by far the most hotspot residues of any gene (120),
reflecting its role as the most frequently mutated tumor suppressor.
```

### Prompt: "Check hotspots for a gene with no known hotspots"

**Tool**: `lookup_cancer_hotspots(gene="FAKEGENE999")`

**Real answer**:
```
Count: 0
(empty list)
```

---

## Driver Mutation Assessment

### Prompt: "Is BRAF V600E a driver mutation?"

**Tool**: `assess_driver_mutation(gene="BRAF", variant="V600E", disease="Melanoma")`

**Expected output** (composite of multiple API calls):
```
Classification: Driver (composite score: 0.85)
Confidence: HIGH

Evidence Signals:
  + Recurrent hotspot (897 samples, q=0.0)          weight: 0.25
  + OncoKB: Oncogenic                                weight: 0.25
  + Strong CIViC evidence (25+ items)                weight: 0.15
  + Oncogenicity SOP: Likely Oncogenic (7+ pts)      weight: 0.20
  + In-silico: Damaging consensus                    weight: 0.08
  + BRAF is a known oncogene                         weight: 0.05
  - gnomAD AF: extremely rare (no penalty)           weight: 0.00

Classification Scale:
  >= 0.70  Driver          — Strong multi-source oncogenic evidence
  0.40-0.69 Likely Driver  — Moderate evidence from 2+ sources
  0.20-0.39 VUS            — Insufficient or conflicting evidence
  < 0.20   Passenger       — Evidence suggests benign/non-functional
```

---

## Cross-Source Observations

### GRCh37 vs GRCh38 Build Sensitivity

| Source | Build-Sensitive? | Notes |
|--------|-----------------|-------|
| gnomAD | Yes | `7-140753336-A-T` (GRCh38) vs `7-140453136-A-T` (GRCh37) — ~300kb offset |
| MyVariant.info | Yes | hgvs_id uses GRCh37 coordinates by default |
| LitVar2 | Returns GRCh38 | `chromosome_position: 7:140753336` from dbSNP |
| Cancer Hotspots | No | Uses protein residue positions (V600), build-agnostic |
| CIViC | No | Queries by gene + variant name |
| ClinVar | No | Queries by gene/variant name or rsID |
| OncoKB | No | Protein-level queries (gene + variant name) |
| UniProt | No | Protein domain positions |
| PubMed | No | Text-based literature search |
| MetaKB | No | Gene + variant name queries |

### Cancer Hotspots API Notes

- **Method**: POST (not GET) — `POST /api/hotspots/single/byGene` with JSON array body `["BRAF"]`
- **q-values**: Returned as strings (`"3.65e-82"`, `"0"`) or `null` — must convert to float
- **3D endpoint**: `POST /api/hotspots/3d/byGene` returns HTTP 500 for most genes (server-side issue, not our bug) — the client handles this gracefully with an empty list fallback
- **Rate limiting**: No documented rate limit, but we apply 5 req/sec as a courtesy

### LitVar2 API Notes

- **Autocomplete**: `GET /variant/autocomplete/?query=BRAF+V600E` — works reliably
- **Variant detail**: `GET /variant/get/{litvar_id}` — **no trailing slash** (returns 404 with trailing slash)
- **`##` in IDs**: LitVar IDs contain `##` suffix (e.g., `litvar@rs113488022##`) — must URL-encode as `%23%23` because `#` is the URL fragment separator
- **v1 Entity**: `GET /api/v1/entity/litvar/rs{id}%23%23` — richer data (diseases, years, all HGVS) from the original LitVar v1 API
- **Publication counts**: LitVar2 autocomplete often reports higher counts than v1 entity (31,630 vs 20,572 for BRAF V600E) — we use the higher count
