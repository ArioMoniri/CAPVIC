"""Unified evidence data models."""

from __future__ import annotations

from datetime import UTC, datetime

from pydantic import BaseModel, Field


class SourceInfo(BaseModel):
    """Information about a data source citation."""

    source_type: str | None = None
    citation: str | None = None
    source_url: str | None = None
    pmid: str | None = None


class CIViCEvidenceItem(BaseModel):
    """A single CIViC evidence item."""

    id: int
    name: str | None = None
    status: str | None = None
    evidence_type: str | None = None
    evidence_level: str | None = None
    evidence_direction: str | None = None
    evidence_rating: int | None = None
    significance: str | None = None
    description: str | None = None
    therapy_interaction_type: str | None = None
    gene: str | None = None
    variant: str | None = None
    molecular_profile: str | None = None
    disease: str | None = None
    disease_doid: str | None = None
    therapies: list[str] = Field(default_factory=list)
    phenotypes: list[str] = Field(default_factory=list)
    source: SourceInfo | None = None


class CIViCAssertion(BaseModel):
    """A CIViC curated assertion."""

    id: int
    name: str | None = None
    assertion_type: str | None = None
    assertion_direction: str | None = None
    significance: str | None = None
    summary: str | None = None
    description: str | None = None
    amp_level: str | None = None
    gene: str | None = None
    variant: str | None = None
    molecular_profile: str | None = None
    disease: str | None = None
    therapies: list[str] = Field(default_factory=list)
    nccn_guideline: str | None = None
    acmg_codes: list[str] = Field(default_factory=list)
    evidence_count: int = 0


class ClinVarVariant(BaseModel):
    """A ClinVar variant record."""

    variation_id: int | None = None
    title: str | None = None
    clinical_significance: str | None = None
    review_status: str | None = None
    review_stars: str | None = None
    genes: list[str] = Field(default_factory=list)
    variation_type: str | None = None
    protein_change: str | None = None
    hgvs_expressions: list[str] = Field(default_factory=list)
    conditions: list[str] = Field(default_factory=list)
    rsid: str | None = None
    last_evaluated: str | None = None
    submitter_count: int = 0
    conflicting: bool = False
    submitter_classifications: list[dict] = Field(default_factory=list)
    clinvar_url: str | None = None


class OncoKBAnnotation(BaseModel):
    """An OncoKB annotation result."""

    oncogenic: str | None = None
    mutation_effect: str | None = None
    known_effect: str | None = None
    highest_sensitive_level: str | None = None
    highest_resistance_level: str | None = None
    highest_diagnostic_level: str | None = None
    highest_prognostic_level: str | None = None
    treatments: list[dict] = Field(default_factory=list)
    gene_summary: str | None = None
    variant_summary: str | None = None
    oncokb_url: str | None = None


class MetaKBInterpretation(BaseModel):
    """A VICC MetaKB interpretation."""

    source: str | None = None
    disease: str | None = None
    drugs: list[str] = Field(default_factory=list)
    evidence_level: str | None = None
    clinical_significance: str | None = None
    description: str | None = None
    url: str | None = None


class GnomADFrequency(BaseModel):
    """gnomAD population allele frequency data."""

    variant_id: str | None = None
    rsid: str | None = None
    allele_frequency: float | None = None
    allele_count: int | None = None
    allele_number: int | None = None
    homozygote_count: int | None = None
    population_frequencies: dict[str, float] = Field(default_factory=dict)
    genome_version: str = "GRCh38"
    filter_status: str | None = None
    source: str = "gnomAD v4"


class Publication(BaseModel):
    """A PubMed publication."""

    pmid: str
    title: str | None = None
    authors: list[str] = Field(default_factory=list)
    journal: str | None = None
    year: str | None = None
    abstract: str | None = None
    mesh_terms: list[str] = Field(default_factory=list)
    doi: str | None = None
    pubmed_url: str | None = None


class PubMedSearchResult(BaseModel):
    """PubMed search result with total count and publications."""

    query: str
    total_count: int = 0
    publications: list[Publication] = Field(default_factory=list)


class EvidenceBundle(BaseModel):
    """Aggregated evidence from all sources for a variant."""

    gene: str
    variant: str
    disease: str | None = None
    retrieved_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    civic_evidence: list[CIViCEvidenceItem] = Field(default_factory=list)
    civic_assertions: list[CIViCAssertion] = Field(default_factory=list)
    clinvar_variants: list[ClinVarVariant] = Field(default_factory=list)
    oncokb_annotation: OncoKBAnnotation | None = None
    metakb_interpretations: list[MetaKBInterpretation] = Field(default_factory=list)
    # Phase 10: scientific enhancement data
    protein_domains: list[ProteinDomain] = Field(default_factory=list)
    in_silico_predictions: InSilicoPredictions | None = None
    gnomad_frequency: GnomADFrequency | None = None
    errors: dict[str, str] = Field(default_factory=dict)

    @property
    def has_civic_data(self) -> bool:
        return bool(self.civic_evidence or self.civic_assertions)

    @property
    def has_clinvar_data(self) -> bool:
        return bool(self.clinvar_variants)

    @property
    def has_oncokb_data(self) -> bool:
        return self.oncokb_annotation is not None

    @property
    def has_gnomad_data(self) -> bool:
        return self.gnomad_frequency is not None

    @property
    def has_domain_data(self) -> bool:
        return bool(self.protein_domains)

    @property
    def has_prediction_data(self) -> bool:
        return self.in_silico_predictions is not None

    @property
    def sources_queried(self) -> list[str]:
        sources = []
        if self.has_civic_data:
            sources.append("CIViC")
        if self.has_clinvar_data:
            sources.append("ClinVar")
        if self.has_oncokb_data:
            sources.append("OncoKB")
        if self.metakb_interpretations:
            sources.append("MetaKB")
        if self.has_gnomad_data:
            sources.append("gnomAD")
        if self.has_domain_data:
            sources.append("UniProt")
        if self.has_prediction_data:
            sources.append("MyVariant.info")
        return sources


# ---------------------------------------------------------------------------
# Variant Normalizer models
# ---------------------------------------------------------------------------


class VariantNotation(BaseModel):
    """Parsed variant notation."""

    original: str
    protein_1letter: str | None = None  # V600E
    protein_3letter: str | None = None  # p.Val600Glu
    cdna: str | None = None  # c.1799T>A
    variant_type: str | None = None  # missense, nonsense, frameshift, splice, etc.
    position: int | None = None  # 600
    ref_aa: str | None = None  # V
    alt_aa: str | None = None  # E


# ---------------------------------------------------------------------------
# UniProt / Protein domain models
# ---------------------------------------------------------------------------


class ProteinDomain(BaseModel):
    """A protein functional domain."""

    name: str
    domain_type: str | None = None  # domain, region, active_site, binding, etc.
    start_pos: int
    end_pos: int
    source: str | None = None  # Pfam, InterPro, PROSITE, etc.
    description: str | None = None


class ProteinFeatures(BaseModel):
    """Protein features from UniProt."""

    gene: str
    uniprot_id: str | None = None
    protein_length: int | None = None
    domains: list[ProteinDomain] = Field(default_factory=list)


class DomainCheckResult(BaseModel):
    """Result of checking whether a variant falls in a functional domain."""

    gene: str
    variant: str
    position: int | None = None
    in_domain: bool = False
    domains: list[ProteinDomain] = Field(default_factory=list)
    evidence_for_om1: bool = False  # True if this supports OM1 evidence code


# ---------------------------------------------------------------------------
# In-silico prediction models
# ---------------------------------------------------------------------------


class InSilicoPredictions(BaseModel):
    """Aggregated in-silico pathogenicity predictions."""

    sift_score: float | None = None
    sift_prediction: str | None = None  # tolerated / deleterious
    polyphen2_score: float | None = None
    polyphen2_prediction: str | None = None  # benign / possibly_damaging / probably_damaging
    revel_score: float | None = None  # meta-predictor, >0.5 likely pathogenic
    cadd_phred: float | None = None  # >20 = top 1% most deleterious
    alphamissense_score: float | None = None
    alphamissense_prediction: str | None = None  # likely_benign / ambiguous / likely_pathogenic
    gerp_score: float | None = None  # >2 = conserved
    phylop_score: float | None = None
    consensus: str | None = None  # "Damaging" / "Benign" / "Mixed"
    damaging_count: int = 0  # How many predictors say damaging
    benign_count: int = 0  # How many say benign
    total_predictors: int = 0

    @property
    def revel_acmg_strength(self) -> str | None:
        """ClinGen SVI-calibrated PP3/BP4 evidence strength from REVEL score.

        Reference: Pejaver et al. 2022, Am J Hum Genet. PMID: 36413997
        Thresholds calibrated by ClinGen Sequence Variant Interpretation WG.
        """
        if self.revel_score is None:
            return None
        if self.revel_score >= 0.773:
            return "PP3_strong"
        if self.revel_score >= 0.644:
            return "PP3_moderate"
        if self.revel_score >= 0.5:
            return "PP3_supporting"
        if self.revel_score <= 0.183:
            return "BP4_strong"
        if self.revel_score <= 0.290:
            return "BP4_moderate"
        if self.revel_score <= 0.4:
            return "BP4_supporting"
        return None  # 0.4 < score < 0.5 — indeterminate
