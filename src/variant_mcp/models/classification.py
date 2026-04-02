"""Classification result models."""

from __future__ import annotations

from pydantic import BaseModel, Field


class AMPTierResult(BaseModel):
    """AMP/ASCO/CAP 4-tier classification result."""

    tier: str = Field(..., description="Tier I, II, III, or IV")
    tier_name: str = Field(..., description="Full tier name")
    evidence_level: str | None = Field(None, description="Level A, B, C, or D (Tier I/II only)")
    confidence: str = Field("LOW", description="LOW, MODERATE, or HIGH")
    evidence_trail: list[str] = Field(default_factory=list)
    sources_used: list[str] = Field(default_factory=list)
    guideline_reference: str = "Li et al., J Mol Diagn 2017. PMID: 27993330"


class OncogenicityScoringResult(BaseModel):
    """ClinGen/CGC/VICC oncogenicity SOP scoring result."""

    total_points: int = 0
    classification: str = "VUS"
    applied_codes: list[AppliedEvidenceCode] = Field(default_factory=list)
    auto_detected: bool = False
    confidence: str = Field("LOW", description="LOW, MODERATE, or HIGH")
    guideline_reference: str = "Horak et al., Genet Med 2022. PMID: 35101336"


class AppliedEvidenceCode(BaseModel):
    """A single applied oncogenicity evidence code."""

    code: str
    points: int
    description: str
    evidence: str | None = None


class ACMGInterpretation(BaseModel):
    """ACMG/AMP germline interpretation summary."""

    aggregate_classification: str | None = None
    review_status: str | None = None
    review_stars: str | None = None
    submitter_count: int = 0
    conflicting_interpretations: bool = False
    submitter_breakdown: dict[str, int] = Field(default_factory=dict)
    last_evaluated: str | None = None
    explanation: str = ""
    criteria_applied: list[str] = Field(default_factory=list)
    guideline_reference: str = "Richards et al., Genet Med 2015. PMID: 25741868"


class VariantClassificationReport(BaseModel):
    """Combined classification report for a variant."""

    gene: str
    variant: str
    disease: str | None = None
    variant_origin: str
    amp_tier: AMPTierResult | None = None
    oncogenicity: OncogenicityScoringResult | None = None
    acmg_interpretation: ACMGInterpretation | None = None
    evidence_summary: dict = Field(default_factory=dict)
    sources_queried: list[str] = Field(default_factory=list)
