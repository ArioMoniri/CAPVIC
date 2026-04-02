"""Pydantic v2 input models for all MCP tools."""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field

from variant_mcp.constants import EvidenceType, SourceDatabase, VariantOrigin


class SearchEvidenceInput(BaseModel):
    """Input for multi-source evidence search."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str = Field(
        ..., min_length=1, max_length=20, description="Gene symbol (e.g., BRAF, KRAS)"
    )
    variant: str | None = Field(
        None, max_length=100, description="Variant name (e.g., V600E, G12C)"
    )
    disease: str | None = Field(None, max_length=200, description="Disease or cancer type")
    therapy: str | None = Field(None, max_length=200, description="Therapy or drug name")
    evidence_type: EvidenceType | None = Field(None, description="Filter by evidence type")
    sources: list[SourceDatabase] | None = Field(
        None, description="Databases to query (default: all available)"
    )
    limit: int = Field(25, ge=1, le=100, description="Max results per source")


class ClassifyVariantInput(BaseModel):
    """Input for variant classification."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str = Field(..., min_length=1, max_length=20, description="Gene symbol")
    variant: str = Field(..., min_length=1, max_length=100, description="Variant name")
    disease: str | None = Field(None, max_length=200, description="Disease or cancer type context")
    variant_origin: VariantOrigin = Field(..., description="somatic or germline")


class GetClinVarInput(BaseModel):
    """Input for ClinVar lookup."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str | None = Field(None, max_length=20, description="Gene symbol")
    variant: str | None = Field(None, max_length=100, description="Variant name")
    variation_id: int | None = Field(None, ge=1, description="ClinVar variation ID")
    rsid: str | None = Field(None, max_length=20, description="dbSNP rsID (e.g., rs113488022)")
    hgvs: str | None = Field(None, max_length=200, description="HGVS expression")
    clinical_significance: str | None = Field(None, description="Filter by significance")


class SearchCIViCInput(BaseModel):
    """Input for CIViC-specific evidence search."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str | None = Field(None, max_length=20, description="Gene symbol")
    variant: str | None = Field(None, max_length=100, description="Variant name")
    disease: str | None = Field(None, max_length=200, description="Disease name")
    therapy: str | None = Field(None, max_length=200, description="Therapy name")
    evidence_type: EvidenceType | None = Field(None, description="Evidence type filter")
    significance: str | None = Field(None, description="Significance filter")
    limit: int = Field(25, ge=1, le=100, description="Max results")


class AnnotateOncoKBInput(BaseModel):
    """Input for OncoKB annotation."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str = Field(..., min_length=1, max_length=20, description="Gene symbol")
    variant: str = Field(..., min_length=1, max_length=100, description="Variant (protein change)")
    tumor_type: str | None = Field(None, max_length=100, description="OncoTree tumor type code")


class SearchMetaKBInput(BaseModel):
    """Input for VICC MetaKB search."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str | None = Field(None, max_length=20, description="Gene symbol")
    variant: str | None = Field(None, max_length=100, description="Variant name")
    disease: str | None = Field(None, max_length=200, description="Disease name")


class LookupGeneInput(BaseModel):
    """Input for gene typeahead lookup."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    query: str = Field(..., min_length=1, max_length=50, description="Gene name or partial name")


class LookupDiseaseInput(BaseModel):
    """Input for disease typeahead lookup."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    query: str = Field(
        ..., min_length=1, max_length=100, description="Disease name or partial name"
    )


class LookupTherapyInput(BaseModel):
    """Input for therapy typeahead lookup."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    query: str = Field(
        ..., min_length=1, max_length=100, description="Therapy/drug name or partial"
    )


class AMPTierInput(BaseModel):
    """Input for AMP/ASCO/CAP tier assignment."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str = Field(..., min_length=1, max_length=20, description="Gene symbol")
    variant: str = Field(..., min_length=1, max_length=100, description="Variant name")
    disease: str | None = Field(None, max_length=200, description="Disease context")


class OncogenicityScoringInput(BaseModel):
    """Input for ClinGen/CGC/VICC oncogenicity SOP scoring."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str = Field(..., min_length=1, max_length=20, description="Gene symbol")
    variant: str = Field(..., min_length=1, max_length=100, description="Variant name")
    evidence_codes: list[str] | None = Field(
        None,
        description="Explicit evidence codes (e.g., ['OVS1', 'OS3']). Auto-detected if omitted.",
    )


class ExplainACMGInput(BaseModel):
    """Input for ACMG criteria explanation."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    criteria_code: str | None = Field(
        None, max_length=10, description="ACMG criteria code (e.g., PVS1, PM2)"
    )
    query: str | None = Field(None, max_length=200, description="General query about ACMG criteria")


class CompareSourcesInput(BaseModel):
    """Input for cross-database variant comparison."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str = Field(..., min_length=1, max_length=20, description="Gene symbol")
    variant: str = Field(..., min_length=1, max_length=100, description="Variant name")


class PathogenicitySummaryInput(BaseModel):
    """Input for pathogenicity summary report."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    gene: str = Field(..., min_length=1, max_length=20, description="Gene symbol")
    variant: str = Field(..., min_length=1, max_length=100, description="Variant name")
    disease: str | None = Field(None, max_length=200, description="Disease context")


class CIViCGetByIdInput(BaseModel):
    """Input for CIViC ID-based lookups."""

    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    id: int = Field(..., ge=1, description="CIViC entity ID")
