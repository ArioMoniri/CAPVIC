"""Tests for Pydantic input and data models."""

import pytest
from pydantic import ValidationError

from variant_mcp.models.classification import (
    AMPTierResult,
    AppliedEvidenceCode,
    OncogenicityScoringResult,
)
from variant_mcp.models.evidence import (
    CIViCEvidenceItem,
    ClinVarVariant,
    EvidenceBundle,
    OncoKBAnnotation,
)
from variant_mcp.models.inputs import (
    AnnotateOncoKBInput,
    ClassifyVariantInput,
    GetClinVarInput,
    OncogenicityScoringInput,
    SearchEvidenceInput,
)


class TestSearchEvidenceInput:
    def test_valid_minimal(self):
        inp = SearchEvidenceInput(gene="BRAF")
        assert inp.gene == "BRAF"
        assert inp.limit == 25

    def test_valid_full(self):
        inp = SearchEvidenceInput(
            gene="KRAS",
            variant="G12C",
            disease="NSCLC",
            therapy="Sotorasib",
            limit=50,
        )
        assert inp.variant == "G12C"

    def test_gene_required(self):
        with pytest.raises(ValidationError):
            SearchEvidenceInput()

    def test_gene_strip_whitespace(self):
        inp = SearchEvidenceInput(gene="  BRAF  ")
        assert inp.gene == "BRAF"

    def test_extra_fields_forbidden(self):
        with pytest.raises(ValidationError):
            SearchEvidenceInput(gene="BRAF", unknown_field="x")

    def test_limit_bounds(self):
        with pytest.raises(ValidationError):
            SearchEvidenceInput(gene="BRAF", limit=0)
        with pytest.raises(ValidationError):
            SearchEvidenceInput(gene="BRAF", limit=101)


class TestClassifyVariantInput:
    def test_valid(self):
        inp = ClassifyVariantInput(gene="BRAF", variant="V600E", variant_origin="somatic")
        assert inp.variant_origin == "somatic"

    def test_germline(self):
        inp = ClassifyVariantInput(gene="BRCA1", variant="R1699W", variant_origin="germline")
        assert inp.variant_origin == "germline"

    def test_invalid_origin(self):
        with pytest.raises(ValidationError):
            ClassifyVariantInput(gene="BRAF", variant="V600E", variant_origin="unknown")


class TestGetClinVarInput:
    def test_variation_id(self):
        inp = GetClinVarInput(variation_id=13961)
        assert inp.variation_id == 13961

    def test_empty_ok(self):
        # All optional — but the client will validate at runtime
        inp = GetClinVarInput()
        assert inp.gene is None


class TestAnnotateOncoKBInput:
    def test_valid(self):
        inp = AnnotateOncoKBInput(gene="BRAF", variant="V600E", tumor_type="MEL")
        assert inp.tumor_type == "MEL"


class TestOncogenicityScoringInput:
    def test_with_codes(self):
        inp = OncogenicityScoringInput(gene="TP53", variant="R175H", evidence_codes=["OVS1", "OS3"])
        assert len(inp.evidence_codes) == 2

    def test_without_codes(self):
        inp = OncogenicityScoringInput(gene="TP53", variant="R175H")
        assert inp.evidence_codes is None


class TestEvidenceBundle:
    def test_empty_bundle(self):
        bundle = EvidenceBundle(gene="BRAF", variant="V600E")
        assert not bundle.has_civic_data
        assert not bundle.has_clinvar_data
        assert not bundle.has_oncokb_data
        assert bundle.sources_queried == []

    def test_with_civic(self):
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            civic_evidence=[
                CIViCEvidenceItem(id=1, evidence_type="PREDICTIVE", evidence_level="A"),
            ],
        )
        assert bundle.has_civic_data
        assert "CIViC" in bundle.sources_queried

    def test_with_clinvar(self):
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            clinvar_variants=[ClinVarVariant(variation_id=13961)],
        )
        assert bundle.has_clinvar_data
        assert "ClinVar" in bundle.sources_queried

    def test_with_oncokb(self):
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            oncokb_annotation=OncoKBAnnotation(oncogenic="Oncogenic"),
        )
        assert bundle.has_oncokb_data
        assert "OncoKB" in bundle.sources_queried


class TestClassificationResults:
    def test_amp_tier_result(self):
        result = AMPTierResult(
            tier="Tier I",
            tier_name="Tier I — Strong Clinical Significance",
            evidence_level="Level A",
            confidence="HIGH",
        )
        assert result.tier == "Tier I"

    def test_oncogenicity_result(self):
        result = OncogenicityScoringResult(
            total_points=14,
            classification="Oncogenic",
            applied_codes=[
                AppliedEvidenceCode(code="OS1", points=4, description="Test"),
                AppliedEvidenceCode(code="OS2", points=4, description="Test"),
                AppliedEvidenceCode(code="OS3", points=4, description="Test"),
                AppliedEvidenceCode(code="OM4", points=2, description="Test"),
            ],
        )
        assert result.classification == "Oncogenic"
        assert result.total_points == 14
