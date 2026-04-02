"""Tests for Phase 10: Scientific Enhancement features."""

from __future__ import annotations

import pytest

from variant_mcp.models.evidence import (
    DomainCheckResult,
    GnomADFrequency,
    InSilicoPredictions,
    ProteinDomain,
    ProteinFeatures,
    Publication,
    PubMedSearchResult,
    VariantNotation,
)
from variant_mcp.utils.variant_normalizer import VariantNormalizer

# ============================================================================
# VariantNormalizer tests
# ============================================================================


class TestVariantNormalizer:
    """Tests for the HGVS notation parser."""

    def setup_method(self) -> None:
        self.normalizer = VariantNormalizer()

    def test_simple_missense(self) -> None:
        result = self.normalizer.normalize("V600E")
        assert result.protein_1letter == "V600E"
        assert result.protein_3letter == "p.Val600Glu"
        assert result.variant_type == "missense"
        assert result.position == 600
        assert result.ref_aa == "V"
        assert result.alt_aa == "E"

    def test_with_p_prefix(self) -> None:
        result = self.normalizer.normalize("p.V600E")
        assert result.protein_1letter == "V600E"
        assert result.protein_3letter == "p.Val600Glu"
        assert result.position == 600

    def test_three_letter_input(self) -> None:
        result = self.normalizer.normalize("p.Val600Glu")
        assert result.protein_1letter == "V600E"
        assert result.variant_type == "missense"
        assert result.position == 600

    def test_nonsense_star(self) -> None:
        result = self.normalizer.normalize("R175*")
        assert result.variant_type == "nonsense"
        assert result.position == 175
        assert result.ref_aa == "R"

    def test_frameshift(self) -> None:
        result = self.normalizer.normalize("K550fs")
        assert result.variant_type == "frameshift"
        assert result.position == 550

    def test_cdna_notation(self) -> None:
        result = self.normalizer.normalize("c.1799T>A")
        assert result.cdna == "c.1799T>A"
        assert result.variant_type == "substitution"

    def test_splice_notation(self) -> None:
        result = self.normalizer.normalize("c.1234+1G>A")
        assert result.variant_type == "splice"

    def test_extract_position(self) -> None:
        assert self.normalizer.extract_position("V600E") == 600
        assert self.normalizer.extract_position("R175*") == 175
        assert self.normalizer.extract_position("p.Val600Glu") == 600

    def test_detect_variant_type_missense(self) -> None:
        assert self.normalizer.detect_variant_type("V600E") == "missense"

    def test_detect_variant_type_nonsense(self) -> None:
        assert self.normalizer.detect_variant_type("R175*") == "nonsense"

    def test_detect_variant_type_frameshift(self) -> None:
        assert self.normalizer.detect_variant_type("K550fs") == "frameshift"

    def test_to_protein_1letter(self) -> None:
        assert self.normalizer.to_protein_1letter("p.Val600Glu") == "V600E"

    def test_to_protein_3letter(self) -> None:
        assert self.normalizer.to_protein_3letter("V600E") == "p.Val600Glu"


# ============================================================================
# New model tests
# ============================================================================


class TestGnomADFrequency:
    def test_default_values(self) -> None:
        freq = GnomADFrequency()
        assert freq.allele_frequency is None
        assert freq.genome_version == "GRCh38"
        assert freq.source == "gnomAD v4"
        assert freq.population_frequencies == {}

    def test_with_data(self) -> None:
        freq = GnomADFrequency(
            variant_id="7-140453136-A-T",
            rsid="rs113488022",
            allele_frequency=0.00001,
            allele_count=3,
            allele_number=300000,
            homozygote_count=0,
            population_frequencies={"nfe": 0.00002, "afr": 0.0},
        )
        assert freq.allele_frequency == 0.00001
        assert freq.population_frequencies["nfe"] == 0.00002


class TestPublication:
    def test_minimal(self) -> None:
        pub = Publication(pmid="27993330")
        assert pub.pmid == "27993330"
        assert pub.authors == []
        assert pub.mesh_terms == []

    def test_with_data(self) -> None:
        pub = Publication(
            pmid="27993330",
            title="Standards and Guidelines",
            authors=["Li MM", "Datto M"],
            journal="J Mol Diagn",
            year="2017",
        )
        assert pub.title == "Standards and Guidelines"
        assert len(pub.authors) == 2


class TestPubMedSearchResult:
    def test_empty(self) -> None:
        result = PubMedSearchResult(query="BRAF V600E")
        assert result.total_count == 0
        assert result.publications == []

    def test_with_results(self) -> None:
        result = PubMedSearchResult(
            query="BRAF V600E",
            total_count=500,
            publications=[Publication(pmid="12345")],
        )
        assert result.total_count == 500
        assert len(result.publications) == 1


class TestVariantNotation:
    def test_parsed_fields(self) -> None:
        notation = VariantNotation(
            original="V600E",
            protein_1letter="V600E",
            protein_3letter="p.Val600Glu",
            variant_type="missense",
            position=600,
            ref_aa="V",
            alt_aa="E",
        )
        assert notation.position == 600
        assert notation.protein_3letter == "p.Val600Glu"


class TestProteinDomain:
    def test_domain(self) -> None:
        domain = ProteinDomain(
            name="Protein kinase domain",
            domain_type="domain",
            start_pos=457,
            end_pos=717,
            source="Pfam",
        )
        assert domain.start_pos == 457
        assert domain.end_pos == 717


class TestProteinFeatures:
    def test_with_domains(self) -> None:
        features = ProteinFeatures(
            gene="BRAF",
            uniprot_id="P15056",
            protein_length=766,
            domains=[
                ProteinDomain(name="Kinase", domain_type="domain", start_pos=457, end_pos=717)
            ],
        )
        assert features.protein_length == 766
        assert len(features.domains) == 1


class TestDomainCheckResult:
    def test_in_domain(self) -> None:
        result = DomainCheckResult(
            gene="BRAF",
            variant="V600E",
            position=600,
            in_domain=True,
            domains=[
                ProteinDomain(name="Kinase", domain_type="domain", start_pos=457, end_pos=717)
            ],
            evidence_for_om1=True,
        )
        assert result.in_domain is True
        assert result.evidence_for_om1 is True

    def test_not_in_domain(self) -> None:
        result = DomainCheckResult(gene="BRAF", variant="A100T", position=100)
        assert result.in_domain is False
        assert result.evidence_for_om1 is False


class TestInSilicoPredictions:
    def test_defaults(self) -> None:
        preds = InSilicoPredictions()
        assert preds.consensus is None
        assert preds.damaging_count == 0
        assert preds.total_predictors == 0

    def test_damaging_consensus(self) -> None:
        preds = InSilicoPredictions(
            sift_score=0.001,
            sift_prediction="deleterious",
            polyphen2_score=0.999,
            polyphen2_prediction="probably_damaging",
            revel_score=0.95,
            cadd_phred=35.0,
            consensus="Damaging",
            damaging_count=4,
            benign_count=0,
            total_predictors=4,
        )
        assert preds.consensus == "Damaging"
        assert preds.damaging_count == 4

    def test_benign_consensus(self) -> None:
        preds = InSilicoPredictions(
            sift_score=0.5,
            sift_prediction="tolerated",
            polyphen2_score=0.1,
            polyphen2_prediction="benign",
            revel_score=0.1,
            consensus="Benign",
            damaging_count=0,
            benign_count=3,
            total_predictors=3,
        )
        assert preds.consensus == "Benign"

    def test_revel_acmg_strength_pp3_strong(self) -> None:
        """REVEL >= 0.773 should yield PP3_strong per ClinGen SVI."""
        preds = InSilicoPredictions(revel_score=0.95)
        assert preds.revel_acmg_strength == "PP3_strong"

    def test_revel_acmg_strength_pp3_moderate(self) -> None:
        """REVEL >= 0.644 should yield PP3_moderate."""
        preds = InSilicoPredictions(revel_score=0.7)
        assert preds.revel_acmg_strength == "PP3_moderate"

    def test_revel_acmg_strength_pp3_supporting(self) -> None:
        """REVEL >= 0.5 should yield PP3_supporting."""
        preds = InSilicoPredictions(revel_score=0.55)
        assert preds.revel_acmg_strength == "PP3_supporting"

    def test_revel_acmg_strength_bp4_strong(self) -> None:
        """REVEL <= 0.183 should yield BP4_strong."""
        preds = InSilicoPredictions(revel_score=0.1)
        assert preds.revel_acmg_strength == "BP4_strong"

    def test_revel_acmg_strength_bp4_moderate(self) -> None:
        """REVEL <= 0.290 should yield BP4_moderate."""
        preds = InSilicoPredictions(revel_score=0.25)
        assert preds.revel_acmg_strength == "BP4_moderate"

    def test_revel_acmg_strength_bp4_supporting(self) -> None:
        """REVEL <= 0.4 should yield BP4_supporting."""
        preds = InSilicoPredictions(revel_score=0.35)
        assert preds.revel_acmg_strength == "BP4_supporting"

    def test_revel_acmg_strength_indeterminate(self) -> None:
        """REVEL between 0.4 and 0.5 should be indeterminate."""
        preds = InSilicoPredictions(revel_score=0.45)
        assert preds.revel_acmg_strength is None

    def test_revel_acmg_strength_none(self) -> None:
        """No REVEL score should return None."""
        preds = InSilicoPredictions()
        assert preds.revel_acmg_strength is None


# ============================================================================
# New input model tests
# ============================================================================


class TestNewInputModels:
    def test_gnomad_input_by_variant_id(self) -> None:
        from variant_mcp.models.inputs import GnomADFrequencyInput

        inp = GnomADFrequencyInput(variant_id="7-140453136-A-T")
        assert inp.variant_id == "7-140453136-A-T"
        assert inp.genome_version == "GRCh38"

    def test_gnomad_input_by_gene_variant(self) -> None:
        from variant_mcp.models.inputs import GnomADFrequencyInput

        inp = GnomADFrequencyInput(gene="BRAF", variant="V600E")
        assert inp.gene == "BRAF"

    def test_normalize_variant_input(self) -> None:
        from variant_mcp.models.inputs import NormalizeVariantInput

        inp = NormalizeVariantInput(variant="V600E", gene="BRAF")
        assert inp.variant == "V600E"

    def test_protein_domain_input(self) -> None:
        from variant_mcp.models.inputs import ProteinDomainInput

        inp = ProteinDomainInput(gene="BRAF", variant="V600E")
        assert inp.gene == "BRAF"

    def test_search_literature_input(self) -> None:
        from variant_mcp.models.inputs import SearchLiteratureInput

        inp = SearchLiteratureInput(gene="BRAF", variant="V600E", disease="Melanoma", limit=5)
        assert inp.limit == 5

    def test_get_publication_input(self) -> None:
        from variant_mcp.models.inputs import GetPublicationInput

        inp = GetPublicationInput(pmid="27993330")
        assert inp.pmid == "27993330"

    def test_predict_variant_effect_input(self) -> None:
        from variant_mcp.models.inputs import PredictVariantEffectInput

        inp = PredictVariantEffectInput(gene="BRAF", variant="V600E")
        assert inp.gene == "BRAF"

    def test_search_literature_limit_bounds(self) -> None:
        from variant_mcp.models.inputs import SearchLiteratureInput

        with pytest.raises(ValueError):
            SearchLiteratureInput(gene="BRAF", limit=100)  # max is 50
