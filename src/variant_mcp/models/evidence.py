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
        return sources
