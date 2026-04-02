"""Constants, API endpoints, enums, evidence codes, and disclaimers."""

from enum import StrEnum

# === API Endpoints ===
CIVIC_GRAPHQL_URL = "https://civicdb.org/api/graphql"
CIVIC_BASE_URL = "https://civicdb.org"
CLINVAR_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CLINVAR_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CLINVAR_BASE_URL = "https://www.ncbi.nlm.nih.gov/clinvar/variation"
ONCOKB_API_URL = "https://www.oncokb.org/api/v1"
METAKB_SEARCH_URL = "https://search.cancervariants.org"

# === Rate Limits (requests per second) ===
CIVIC_RATE_LIMIT = 3
CLINVAR_RATE_LIMIT = 3
ONCOKB_RATE_LIMIT = 5
METAKB_RATE_LIMIT = 3

# === Timeouts ===
REQUEST_TIMEOUT = 30
MAX_RETRIES = 3
RETRY_BACKOFF_BASE = 1.0


# === Enums ===
class VariantOrigin(StrEnum):
    SOMATIC = "somatic"
    GERMLINE = "germline"


class EvidenceType(StrEnum):
    PREDICTIVE = "PREDICTIVE"
    DIAGNOSTIC = "DIAGNOSTIC"
    PROGNOSTIC = "PROGNOSTIC"
    PREDISPOSING = "PREDISPOSING"
    ONCOGENIC = "ONCOGENIC"
    FUNCTIONAL = "FUNCTIONAL"


class EvidenceSignificance(StrEnum):
    SENSITIVITY = "SENSITIVITY"
    RESISTANCE = "RESISTANCE"
    POSITIVE = "POSITIVE"
    NEGATIVE = "NEGATIVE"
    BETTER_OUTCOME = "BETTER_OUTCOME"
    POOR_OUTCOME = "POOR_OUTCOME"
    PATHOGENIC = "PATHOGENIC"
    LIKELY_PATHOGENIC = "LIKELY_PATHOGENIC"
    BENIGN = "BENIGN"
    LIKELY_BENIGN = "LIKELY_BENIGN"
    UNCERTAIN_SIGNIFICANCE = "UNCERTAIN_SIGNIFICANCE"
    ONCOGENIC = "ONCOGENIC"
    LIKELY_ONCOGENIC = "LIKELY_ONCOGENIC"


class SourceDatabase(StrEnum):
    CIVIC = "civic"
    CLINVAR = "clinvar"
    ONCOKB = "oncokb"
    METAKB = "metakb"


# === AMP/ASCO/CAP 4-Tier Classification (2017 + 2025 update) ===
AMP_TIER_DEFINITIONS = {
    "I": {
        "name": "Tier I — Strong Clinical Significance",
        "level_A": (
            "FDA-approved therapy, included in professional guidelines "
            "(NCCN, CAP, ASCO, ESMO), or well-powered studies with consensus"
        ),
        "level_B": (
            "Well-powered studies with expert consensus from clinical trials or large cohorts"
        ),
    },
    "II": {
        "name": "Tier II — Potential Clinical Significance",
        "level_C": (
            "FDA-approved therapies for different tumor types, or investigational "
            "therapies with evidence from clinical trials"
        ),
        "level_D": (
            "Preclinical trials or case reports showing potential association "
            "with diagnosis, prognosis, or therapy"
        ),
    },
    "III": {
        "name": "Tier III — Unknown Clinical Significance",
        "description": (
            "Variant not observed at a significant frequency in population databases, "
            "no convincing evidence of cancer association published"
        ),
    },
    "IV": {
        "name": "Tier IV — Benign or Likely Benign",
        "description": (
            "Observed at significant frequency in population databases (gnomAD/ExAC). "
            "No published evidence of cancer association."
        ),
    },
}

# === ClinGen/CGC/VICC Oncogenicity SOP Evidence Codes ===
# Reference: Horak et al., Genetics in Medicine 2022
ONCOGENICITY_EVIDENCE_CODES = {
    # Oncogenic codes (positive evidence)
    "OVS1": {
        "strength": "very_strong",
        "points": 8,
        "description": "Null variant in bona fide tumor suppressor gene",
    },
    "OS1": {
        "strength": "strong",
        "points": 4,
        "description": "Same amino acid change as established oncogenic variant",
    },
    "OS2": {
        "strength": "strong",
        "points": 4,
        "description": "Well-established functional studies show oncogenic effect",
    },
    "OS3": {
        "strength": "strong",
        "points": 4,
        "description": "Located in hotspot with sufficient statistical evidence",
    },
    "OM1": {
        "strength": "moderate",
        "points": 2,
        "description": "Located in critical/functional domain without benign variation",
    },
    "OM2": {
        "strength": "moderate",
        "points": 2,
        "description": "Protein length changes in known oncogene/TSG",
    },
    "OM3": {
        "strength": "moderate",
        "points": 2,
        "description": "At same position as established oncogenic missense",
    },
    "OM4": {
        "strength": "moderate",
        "points": 2,
        "description": "Absent/extremely rare in population databases",
    },
    "OP1": {
        "strength": "supporting",
        "points": 1,
        "description": "All computational evidence supports oncogenic effect",
    },
    "OP2": {
        "strength": "supporting",
        "points": 1,
        "description": "Somatic variant in gene with single cancer etiology",
    },
    "OP3": {
        "strength": "supporting",
        "points": 1,
        "description": "Located in hotspot with limited statistical evidence",
    },
    "OP4": {
        "strength": "supporting",
        "points": 1,
        "description": "Absent in population databases",
    },
    # Benign codes (negative evidence)
    "SBVS1": {
        "strength": "very_strong_benign",
        "points": -8,
        "description": "MAF >5% in any general continental population in gnomAD",
    },
    "SBS1": {
        "strength": "strong_benign",
        "points": -4,
        "description": "MAF >1% in any general continental population in gnomAD",
    },
    "SBS2": {
        "strength": "strong_benign",
        "points": -4,
        "description": "Well-established functional studies show no oncogenic effect",
    },
    "SBP1": {
        "strength": "supporting_benign",
        "points": -1,
        "description": "All computational evidence suggests benign",
    },
    "SBP2": {
        "strength": "supporting_benign",
        "points": -1,
        "description": "Synonymous change with no predicted splice impact",
    },
}

# Oncogenicity classification thresholds (point-based)
ONCOGENICITY_THRESHOLDS = {
    "Oncogenic": 10,
    "Likely Oncogenic": 6,
    "VUS": 0,  # between -1 and 5 inclusive
    "Likely Benign": -6,
    "Benign": -10,
}

# === ACMG/AMP Germline Classification Terms ===
ACMG_CLASSIFICATIONS = [
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
]

# ACMG/AMP Evidence Criteria
ACMG_CRITERIA = {
    # Very Strong
    "PVS1": "Null variant in a gene where LOF is a known mechanism of disease",
    # Strong Pathogenic
    "PS1": "Same amino acid change as a previously established pathogenic variant",
    "PS2": "De novo (confirmed) in a patient with disease, no family history",
    "PS3": "Well-established in vitro or in vivo functional studies support damaging effect",
    "PS4": "Prevalence in affected individuals significantly increased vs controls",
    # Moderate Pathogenic
    "PM1": "Located in a mutational hot spot and/or critical functional domain",
    "PM2": "Absent from controls (or at extremely low frequency) in population databases",
    "PM3": "For recessive disorders, detected in trans with a pathogenic variant",
    "PM4": "Protein length changes as a result of in-frame deletions/insertions",
    "PM5": "Novel missense change at an amino acid where a different pathogenic missense has been seen",
    "PM6": "Assumed de novo, but without confirmation of paternity and maternity",
    # Supporting Pathogenic
    "PP1": "Cosegregation with disease in multiple affected family members",
    "PP2": "Missense variant in gene with low rate of benign missense variation",
    "PP3": "Multiple lines of computational evidence support a deleterious effect",
    "PP4": "Patient's phenotype or family history is highly specific for disease with single genetic etiology",
    "PP5": "Reputable source recently reports variant as pathogenic (with evidence)",
    # Benign Standalone
    "BA1": "Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or gnomAD",
    # Strong Benign
    "BS1": "Allele frequency is greater than expected for disorder",
    "BS2": "Observed in a healthy adult individual for a recessive, dominant, or X-linked disorder",
    "BS3": "Well-established in vitro or in vivo functional studies show no damaging effect",
    "BS4": "Lack of segregation in affected members of a family",
    # Supporting Benign
    "BP1": "Missense variant in a gene for which primarily truncating variants are known to cause disease",
    "BP2": "Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder",
    "BP3": "In-frame deletions/insertions in a repetitive region without a known function",
    "BP4": "Multiple lines of computational evidence suggest no impact on gene or gene product",
    "BP5": "Variant found in a case with an alternate molecular basis for disease",
    "BP6": "Reputable source recently reports variant as benign (with evidence)",
    "BP7": "A synonymous variant for which splicing prediction algorithms predict no impact",
}

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

# Known tumor suppressor genes (for OVS1 code auto-detection)
KNOWN_TUMOR_SUPPRESSORS = {
    "TP53",
    "RB1",
    "BRCA1",
    "BRCA2",
    "APC",
    "PTEN",
    "VHL",
    "WT1",
    "NF1",
    "NF2",
    "CDKN2A",
    "SMAD4",
    "STK11",
    "MLH1",
    "MSH2",
    "MSH6",
    "PMS2",
    "CDH1",
    "PTCH1",
    "MEN1",
    "TSC1",
    "TSC2",
    "BMPR1A",
    "SMARCB1",
    "BAP1",
    "ARID1A",
    "KDM6A",
    "FBXW7",
    "RNF43",
    "ZNRF3",
    "LATS1",
    "LATS2",
    "FAT1",
    "AXIN1",
    "AXIN2",
}

# Known oncogenes (for classification context)
KNOWN_ONCOGENES = {
    "BRAF",
    "KRAS",
    "NRAS",
    "HRAS",
    "EGFR",
    "PIK3CA",
    "ALK",
    "RET",
    "MET",
    "KIT",
    "PDGFRA",
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "ABL1",
    "JAK2",
    "MPL",
    "ERBB2",
    "ERBB3",
    "IDH1",
    "IDH2",
    "CTNNB1",
    "SMO",
    "GNA11",
    "GNAQ",
    "AKT1",
    "MAP2K1",
    "MAP2K2",
    "MTOR",
    "MYC",
    "MYCN",
    "MYCL",
    "CDK4",
    "CDK6",
    "CCND1",
    "CCND2",
    "CCND3",
    "MDM2",
    "MDM4",
}

# Null variant types (for OVS1 auto-detection)
NULL_VARIANT_TYPES = {
    "nonsense",
    "frameshift",
    "splice_site",
    "splice_donor",
    "splice_acceptor",
    "start_loss",
    "exon_deletion",
    "transcript_ablation",
}

DISCLAIMER = (
    "**RESEARCH DISCLAIMER**: This tool aggregates data from CIViC, ClinVar, "
    "OncoKB, and VICC MetaKB for **research and educational purposes only**. "
    "Classification suggestions follow published AMP/ASCO/CAP, ClinGen/CGC/VICC, "
    "and ACMG/AMP frameworks but are **computationally derived approximations** "
    "— NOT validated clinical interpretations. Clinical variant classification "
    "MUST be performed by qualified molecular pathologists using validated "
    "laboratory procedures. Do not use for clinical decision-making without "
    "professional review."
)

# Key reference PMIDs
REFERENCE_PMIDS = {
    "AMP_ASCO_CAP": "27993330",  # Li et al. 2017
    "ONCOGENICITY_SOP": "35101336",  # Horak et al. 2022
    "ACMG_AMP": "25741868",  # Richards et al. 2015
    "CIVIC": "28138153",  # Griffith et al. 2017
    "VICC_METAKB": "32246132",  # Wagner et al. 2020
}
