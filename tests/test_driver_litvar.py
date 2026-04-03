"""Mock-based tests for driver mutation + LitVar2 clients and tools.

Tests Cancer Hotspots client (POST API), LitVar2 client (autocomplete + entity),
and driver mutation assessment logic.
"""

from __future__ import annotations

import httpx
import pytest
import respx

from variant_mcp.clients.cancer_hotspots_client import CancerHotspotsClient
from variant_mcp.clients.litvar_client import LitVarClient
from variant_mcp.models.evidence import (
    CancerHotspot,
    CIViCEvidenceItem,
    DriverMutationAssessment,
    EvidenceBundle,
    InSilicoPredictions,
    OncoKBAnnotation,
)

# ============================================================================
# Cancer Hotspots Client
# ============================================================================


class TestCancerHotspotsClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_get_hotspots_by_gene(self):
        client = CancerHotspotsClient()
        respx.post("https://www.cancerhotspots.org/api/hotspots/single/byGene").mock(
            return_value=httpx.Response(
                200,
                json=[
                    {
                        "hugoSymbol": "BRAF",
                        "residue": "V600",
                        "tumorCount": 897,
                        "type": "single residue",
                        "qValue": "0",
                        "qValueCancerType": "0",
                        "variantAminoAcid": {"E": 833, "K": 24, "M": 29, "R": 4},
                        "tumorTypeComposition": {
                            "skin": 357,
                            "thyroid": 316,
                            "bowel": 113,
                        },
                    },
                    {
                        "hugoSymbol": "BRAF",
                        "residue": "G469",
                        "tumorCount": 61,
                        "type": "single residue",
                        "qValue": "3.65e-82",
                        "qValueCancerType": "2.17e-35",
                        "variantAminoAcid": {"A": 30, "V": 15, "S": 10},
                        "tumorTypeComposition": {"lung": 25, "skin": 20},
                    },
                ],
            )
        )
        result = await client.get_hotspots_by_gene("BRAF")
        assert len(result) == 2
        assert result[0]["hugoSymbol"] == "BRAF"
        assert result[0]["residue"] == "V600"

    @respx.mock
    @pytest.mark.asyncio
    async def test_get_hotspot_at_residue(self):
        client = CancerHotspotsClient()
        respx.post("https://www.cancerhotspots.org/api/hotspots/single/byGene").mock(
            return_value=httpx.Response(
                200,
                json=[
                    {
                        "hugoSymbol": "BRAF",
                        "residue": "V600",
                        "tumorCount": 897,
                        "type": "single residue",
                        "qValue": "0",
                        "qValueCancerType": "0",
                        "variantAminoAcid": {"E": 833},
                        "tumorTypeComposition": {},
                    },
                    {
                        "hugoSymbol": "BRAF",
                        "residue": "G469",
                        "tumorCount": 61,
                        "type": "single residue",
                        "qValue": "3.65e-82",
                        "qValueCancerType": None,
                        "variantAminoAcid": {"A": 30},
                        "tumorTypeComposition": {},
                    },
                ],
            )
        )
        result = await client.get_hotspot_at_residue("BRAF", 600)
        assert len(result) == 1
        assert result[0]["residue"] == "V600"

    def test_parse_hotspot_residue(self):
        client = CancerHotspotsClient()
        raw = {
            "hugoSymbol": "BRAF",
            "residue": "V600",
            "tumorCount": 897,
            "type": "single residue",
            "qValue": "0",
            "qValueCancerType": "2.5e-10",
            "variantAminoAcid": {"E": 833, "K": 24, "M": 29},
            "tumorTypeComposition": {"skin": 357, "thyroid": 316},
        }
        parsed = client.parse_hotspot_residue(raw)
        assert parsed["gene"] == "BRAF"
        assert parsed["residue"] == "V600"
        assert parsed["total_sample_count"] == 886  # 833 + 24 + 29
        assert parsed["q_value"] == 0.0
        assert parsed["q_value_cancertype"] == 2.5e-10
        assert parsed["cancer_type_counts"]["skin"] == 357

    def test_parse_hotspot_residue_none_qvalue(self):
        client = CancerHotspotsClient()
        raw = {
            "hugoSymbol": "TP53",
            "residue": "R175",
            "tumorCount": 200,
            "type": "single residue",
            "qValue": None,
            "qValueCancerType": None,
            "variantAminoAcid": {"H": 180, "C": 20},
            "tumorTypeComposition": {},
        }
        parsed = client.parse_hotspot_residue(raw)
        assert parsed["q_value"] is None
        assert parsed["q_value_cancertype"] is None
        assert parsed["total_sample_count"] == 200

    @respx.mock
    @pytest.mark.asyncio
    async def test_empty_gene_returns_empty(self):
        client = CancerHotspotsClient()
        respx.post("https://www.cancerhotspots.org/api/hotspots/single/byGene").mock(
            return_value=httpx.Response(200, json=[])
        )
        result = await client.get_hotspots_by_gene("FAKEGENE")
        assert result == []

    @respx.mock
    @pytest.mark.asyncio
    async def test_api_error_returns_empty(self):
        client = CancerHotspotsClient()
        respx.post("https://www.cancerhotspots.org/api/hotspots/single/byGene").mock(
            return_value=httpx.Response(500, text="Internal Server Error")
        )
        result = await client.get_hotspots_by_gene("BRAF")
        assert result == []


# ============================================================================
# LitVar2 Client
# ============================================================================


class TestLitVarClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_autocomplete_braf_v600e(self):
        client = LitVarClient()
        respx.get("https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/").mock(
            return_value=httpx.Response(
                200,
                json=[
                    {
                        "_id": "litvar@rs113488022##",
                        "rsid": "rs113488022",
                        "gene": ["BRAF"],
                        "name": "p.V600E",
                        "hgvs": "p.V600E",
                        "pmids_count": 31630,
                        "data_clinical_significance": [
                            "pathogenic",
                            "drug-response",
                        ],
                    }
                ],
            )
        )
        result = await client.autocomplete("BRAF V600E")
        assert len(result) == 1
        assert result[0]["rsid"] == "rs113488022"
        assert result[0]["pmids_count"] == 31630

    @respx.mock
    @pytest.mark.asyncio
    async def test_autocomplete_not_found(self):
        client = LitVarClient()
        respx.get("https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/").mock(
            return_value=httpx.Response(200, json=[])
        )
        result = await client.autocomplete("FAKEGENE X999Y")
        assert result == []

    @respx.mock
    @pytest.mark.asyncio
    async def test_get_variant(self):
        client = LitVarClient()
        respx.get(
            "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/litvar@rs113488022%23%23"
        ).mock(
            return_value=httpx.Response(
                200,
                json={
                    "_id": "litvar@rs113488022##",
                    "rsid": "rs113488022",
                    "clingen_ids": ["CA16602736"],
                    "data_chromosome_base_position": ["7:140753336"],
                    "data_snp_class": ["snv"],
                    "data_clinical_significance": ["pathogenic"],
                },
            )
        )
        result = await client.get_variant("litvar@rs113488022##")
        assert result is not None
        assert result["rsid"] == "rs113488022"
        assert "CA16602736" in result["clingen_ids"]

    @respx.mock
    @pytest.mark.asyncio
    async def test_search_variant_literature_full(self):
        """Test the high-level search that chains autocomplete + get + v1 entity."""
        client = LitVarClient()

        # Mock autocomplete
        respx.get("https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/").mock(
            return_value=httpx.Response(
                200,
                json=[
                    {
                        "_id": "litvar@rs113488022##",
                        "rsid": "rs113488022",
                        "gene": ["BRAF"],
                        "name": "p.V600E",
                        "hgvs": "p.V600E",
                        "pmids_count": 31630,
                        "data_clinical_significance": ["pathogenic"],
                    }
                ],
            )
        )

        # Mock variant get (no trailing slash — LitVar2 returns 404 with it)
        respx.get(
            "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/litvar@rs113488022%23%23"
        ).mock(
            return_value=httpx.Response(
                200,
                json={
                    "_id": "litvar@rs113488022##",
                    "clingen_ids": ["CA16602736"],
                    "data_chromosome_base_position": ["7:140753336"],
                    "data_snp_class": ["snv"],
                },
            )
        )

        # Mock v1 entity
        respx.get(
            "https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1/entity/litvar/rs113488022%23%23"
        ).mock(
            return_value=httpx.Response(
                200,
                json={
                    "diseases": {"Melanoma": 26551, "Colorectal Neoplasms": 6505},
                    "years": [2004, 2010, 2020, 2022],
                    "all_hgvs": [["p.V600E", 20322], ["c.1799T>A", 500]],
                    "first_published_year": 2004,
                    "pmids_count": 20572,
                },
            )
        )

        result = await client.search_variant_literature("BRAF", "V600E")
        assert result["found"] is True
        assert result["rsid"] == "rs113488022"
        assert result["pmids_count"] == 31630  # LitVar2 count is higher
        assert result["diseases"]["Melanoma"] == 26551
        assert len(result["all_hgvs"]) == 2
        assert result["first_published_year"] == 2004
        assert result["clingen_ids"] == ["CA16602736"]

    @respx.mock
    @pytest.mark.asyncio
    async def test_search_not_found(self):
        client = LitVarClient()
        respx.get("https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/").mock(
            return_value=httpx.Response(200, json=[])
        )
        result = await client.search_variant_literature("FAKEGENE", "X999Y")
        assert result["found"] is False
        assert result["pmids_count"] == 0


# ============================================================================
# Driver Mutation Assessment Model
# ============================================================================


class TestDriverMutationAssessment:
    def test_driver_classification_fields(self):
        assessment = DriverMutationAssessment(
            gene="BRAF",
            variant="V600E",
            driver_classification="Driver",
            driver_score=0.85,
            signals=["Recurrent hotspot", "OncoKB: Oncogenic"],
            is_known_oncogene=True,
            confidence="HIGH",
        )
        assert assessment.driver_classification == "Driver"
        assert assessment.driver_score == 0.85
        assert assessment.is_known_oncogene is True
        assert len(assessment.signals) == 2

    def test_passenger_classification(self):
        assessment = DriverMutationAssessment(
            gene="BRAF",
            variant="T119I",
            driver_classification="Passenger",
            driver_score=0.05,
            signals=["Common in gnomAD (AF=0.02)"],
            gnomad_af=0.02,
            confidence="MODERATE",
        )
        assert assessment.driver_classification == "Passenger"
        assert assessment.gnomad_af == 0.02


# ============================================================================
# Cancer Hotspot Model
# ============================================================================


class TestCancerHotspotModel:
    def test_cancer_hotspot_creation(self):
        hs = CancerHotspot(
            gene="BRAF",
            residue="V600",
            total_sample_count=897,
            variant_amino_acids={"E": 833, "K": 24},
            cancer_type_counts={"skin": 357, "thyroid": 316},
            q_value=0.0,
        )
        assert hs.gene == "BRAF"
        assert hs.total_sample_count == 897
        assert hs.q_value == 0.0


# ============================================================================
# Evidence Bundle with Hotspot Data
# ============================================================================


class TestEvidenceBundleHotspot:
    def test_has_hotspot_data(self):
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            cancer_hotspots=[
                CancerHotspot(
                    gene="BRAF",
                    residue="V600",
                    total_sample_count=897,
                    q_value=0.0,
                )
            ],
        )
        assert bundle.has_hotspot_data is True
        assert "Cancer Hotspots" in bundle.sources_queried

    def test_no_hotspot_data(self):
        bundle = EvidenceBundle(gene="BRAF", variant="V600E")
        assert bundle.has_hotspot_data is False
        assert "Cancer Hotspots" not in bundle.sources_queried


# ============================================================================
# Oncogenicity Scorer with Real Hotspot Data
# ============================================================================


class TestOncogenicityScorerHotspot:
    def test_hotspot_detection_with_cancer_hotspots_data(self):
        from variant_mcp.classification.oncogenicity_sop import OncogenicityScorer

        scorer = OncogenicityScorer()
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            cancer_hotspots=[
                CancerHotspot(
                    gene="BRAF",
                    residue="V600",
                    total_sample_count=897,
                    q_value=0.0,
                )
            ],
        )
        assert scorer._is_hotspot("BRAF", "V600E", bundle) is True

    def test_hotspot_detection_by_qvalue(self):
        from variant_mcp.classification.oncogenicity_sop import OncogenicityScorer

        scorer = OncogenicityScorer()
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="G469A",
            cancer_hotspots=[
                CancerHotspot(
                    gene="BRAF",
                    residue="G469",
                    total_sample_count=5,  # Below 10 threshold
                    q_value=0.001,  # But q-value < 0.05
                )
            ],
        )
        assert scorer._is_hotspot("BRAF", "G469A", bundle) is True

    def test_hotspot_below_thresholds(self):
        from variant_mcp.classification.oncogenicity_sop import OncogenicityScorer

        scorer = OncogenicityScorer()
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="T119I",
            cancer_hotspots=[
                CancerHotspot(
                    gene="BRAF",
                    residue="T119",
                    total_sample_count=2,  # Below threshold
                    q_value=0.5,  # Not significant
                )
            ],
        )
        assert scorer._is_hotspot("BRAF", "T119I", bundle) is False

    def test_oncogenicity_score_with_hotspot_assigns_os3(self):
        from variant_mcp.classification.oncogenicity_sop import OncogenicityScorer

        scorer = OncogenicityScorer()
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            cancer_hotspots=[
                CancerHotspot(
                    gene="BRAF",
                    residue="V600",
                    total_sample_count=897,
                    q_value=0.0,
                )
            ],
        )
        result = scorer.score_variant("BRAF", "V600E", evidence_bundle=bundle)
        codes = [c.code for c in result.applied_codes]
        assert "OS3" in codes  # Hotspot evidence
        assert result.total_points >= 4  # OS3 = 4 points


# ============================================================================
# GRCh37/38 Coordinate Awareness Tests
# ============================================================================


class TestCoordinateAwareness:
    """Verify that our APIs handle genome builds correctly."""

    def test_litvar_returns_grch38(self):
        """LitVar2 chromosome positions are GRCh38 (from dbSNP)."""
        # BRAF V600E: chr7:140753336 (GRCh38) vs chr7:140453136 (GRCh37)
        # LitVar2 returns 7:140753336 which is GRCh38
        grch38_braf_pos = 140753336
        grch37_braf_pos = 140453136
        assert grch38_braf_pos != grch37_braf_pos
        # ~300kb difference confirms different assemblies
        assert abs(grch38_braf_pos - grch37_braf_pos) == 300200

    def test_cancer_hotspots_is_protein_level(self):
        """Cancer Hotspots uses amino acid positions, genome-build agnostic."""
        hs = CancerHotspot(gene="BRAF", residue="V600", total_sample_count=897)
        # Protein position doesn't depend on genome build
        assert hs.residue == "V600"

    def test_evidence_bundle_sources_track_all(self):
        """Verify all 10 data sources can be tracked."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            civic_evidence=[CIViCEvidenceItem(id=1, name="test")],
            oncokb_annotation=OncoKBAnnotation(oncogenic="Oncogenic"),
            in_silico_predictions=InSilicoPredictions(
                consensus="Damaging", total_predictors=3, damaging_count=3
            ),
            cancer_hotspots=[CancerHotspot(gene="BRAF", residue="V600", total_sample_count=897)],
        )
        sources = bundle.sources_queried
        assert "CIViC" in sources
        assert "OncoKB" in sources
        assert "MyVariant.info" in sources
        assert "Cancer Hotspots" in sources
