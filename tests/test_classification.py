"""Tests for classification engines."""

from variant_mcp.classification.acmg_amp import ACMGAMPHelper
from variant_mcp.classification.amp_asco_cap import AMPTierClassifier
from variant_mcp.classification.oncogenicity_sop import OncogenicityScorer
from variant_mcp.models.evidence import (
    CIViCEvidenceItem,
    ClinVarVariant,
    EvidenceBundle,
    InSilicoPredictions,
    OncoKBAnnotation,
    ProteinDomain,
)


class TestAMPTierClassifier:
    def setup_method(self):
        self.classifier = AMPTierClassifier()

    def test_empty_evidence_tier_iii(self):
        """No evidence should result in Tier III."""
        bundle = EvidenceBundle(gene="GENE1", variant="X123Y")
        result = self.classifier.classify(bundle)
        assert result.tier == "Tier III"
        assert result.confidence == "LOW"

    def test_oncokb_level1_tier_ia(self):
        """OncoKB Level 1 should map to Tier I/A."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            oncokb_annotation=OncoKBAnnotation(
                oncogenic="Oncogenic",
                known_effect="Gain-of-function",
                highest_sensitive_level="LEVEL_1",
            ),
        )
        result = self.classifier.classify(bundle)
        assert result.tier == "Tier I"
        assert result.evidence_level == "Level A"

    def test_oncokb_level2_tier_ib(self):
        """OncoKB Level 2 should map to Tier I/B."""
        bundle = EvidenceBundle(
            gene="KRAS",
            variant="G12C",
            oncokb_annotation=OncoKBAnnotation(
                oncogenic="Oncogenic",
                highest_sensitive_level="LEVEL_2",
            ),
        )
        result = self.classifier.classify(bundle)
        assert result.tier == "Tier I"
        assert result.evidence_level == "Level B"

    def test_oncokb_level3a_tier_iic(self):
        """OncoKB Level 3A should map to Tier II/C."""
        bundle = EvidenceBundle(
            gene="PIK3CA",
            variant="H1047R",
            oncokb_annotation=OncoKBAnnotation(
                oncogenic="Oncogenic",
                highest_sensitive_level="LEVEL_3A",
            ),
        )
        result = self.classifier.classify(bundle)
        assert result.tier == "Tier II"
        assert result.evidence_level == "Level C"

    def test_civic_level_a_tier_ia(self):
        """CIViC Level A evidence should push to Tier I/A."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            civic_evidence=[
                CIViCEvidenceItem(
                    id=1,
                    evidence_type="PREDICTIVE",
                    evidence_level="A",
                    significance="SENSITIVITY",
                ),
                CIViCEvidenceItem(
                    id=2,
                    evidence_type="PREDICTIVE",
                    evidence_level="A",
                    significance="SENSITIVITY",
                ),
            ],
        )
        result = self.classifier.classify(bundle)
        assert result.tier == "Tier I"
        assert result.evidence_level == "Level A"

    def test_clinvar_benign_tier_iv(self):
        """ClinVar benign classification should support Tier IV."""
        bundle = EvidenceBundle(
            gene="GENE1",
            variant="X123Y",
            clinvar_variants=[
                ClinVarVariant(
                    variation_id=1,
                    clinical_significance="Benign",
                    review_status="criteria provided, multiple submitters, no conflicts",
                    review_stars="★★",
                ),
            ],
        )
        result = self.classifier.classify(bundle)
        assert result.tier == "Tier IV"

    def test_multi_source_high_confidence(self):
        """Multiple concordant sources should give HIGH confidence."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            oncokb_annotation=OncoKBAnnotation(
                oncogenic="Oncogenic",
                highest_sensitive_level="LEVEL_1",
            ),
            civic_evidence=[
                CIViCEvidenceItem(id=1, evidence_type="PREDICTIVE", evidence_level="A"),
            ],
            clinvar_variants=[
                ClinVarVariant(variation_id=1, clinical_significance="Pathogenic"),
            ],
        )
        result = self.classifier.classify(bundle)
        assert result.tier == "Tier I"
        assert result.confidence == "HIGH"


class TestOncogenicityScorer:
    def setup_method(self):
        self.scorer = OncogenicityScorer()

    def test_explicit_codes_oncogenic(self):
        """Explicit codes summing >= 10 should classify as Oncogenic."""
        result = self.scorer.score_variant("BRAF", "V600E", evidence_codes=["OS1", "OS2", "OS3"])
        assert result.total_points == 12
        assert result.classification == "Oncogenic"
        assert not result.auto_detected

    def test_explicit_codes_likely_oncogenic(self):
        """Codes summing 6-9 should classify as Likely Oncogenic."""
        result = self.scorer.score_variant("BRAF", "V600E", evidence_codes=["OS1", "OM1"])
        assert result.total_points == 6
        assert result.classification == "Likely Oncogenic"

    def test_explicit_codes_vus(self):
        """Codes summing 0-5 should classify as VUS."""
        result = self.scorer.score_variant("GENE1", "X123Y", evidence_codes=["OP1"])
        assert result.total_points == 1
        assert result.classification == "VUS"

    def test_explicit_codes_benign(self):
        """Codes summing <= -10 should classify as Benign."""
        result = self.scorer.score_variant("GENE1", "X123Y", evidence_codes=["SBVS1", "SBS1"])
        assert result.total_points == -12
        assert result.classification == "Benign"

    def test_explicit_codes_likely_benign(self):
        """Codes summing -6 to -9 should classify as Likely Benign."""
        result = self.scorer.score_variant("GENE1", "X123Y", evidence_codes=["SBVS1"])
        assert result.total_points == -8
        assert result.classification == "Likely Benign"

    def test_auto_detect_null_variant_tsg(self):
        """Null variant in TSG should auto-detect OVS1."""
        result = self.scorer.score_variant("TP53", "R175*")
        assert result.auto_detected
        codes = [c.code for c in result.applied_codes]
        assert "OVS1" in codes

    def test_auto_detect_frameshift_tsg(self):
        """Frameshift in TSG should auto-detect OVS1."""
        result = self.scorer.score_variant("BRCA1", "E23fs")
        codes = [c.code for c in result.applied_codes]
        assert "OVS1" in codes

    def test_explain_codes(self):
        """Code explanations should be human-readable."""
        explanation = self.scorer.explain_codes(["OVS1", "OS3", "SBS1"])
        assert "OVS1" in explanation
        assert "OS3" in explanation
        assert "SBS1" in explanation
        assert "points" in explanation

    def test_unknown_code_warning(self):
        """Unknown codes should be silently ignored."""
        result = self.scorer.score_variant("GENE1", "X123Y", evidence_codes=["INVALID_CODE"])
        assert result.total_points == 0

    def test_auto_detect_om1_with_uniprot_domains(self):
        """OM1 should use real UniProt domain data when available."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            protein_domains=[
                ProteinDomain(
                    name="Protein kinase", domain_type="Domain", start_pos=457, end_pos=717
                )
            ],
        )
        result = self.scorer.score_variant("BRAF", "V600E", evidence_bundle=bundle)
        codes = [c.code for c in result.applied_codes]
        assert "OM1" in codes
        om1 = next(c for c in result.applied_codes if c.code == "OM1")
        assert "UniProt" in om1.evidence

    def test_auto_detect_op1_with_damaging_predictions(self):
        """OP1 should fire when in-silico consensus is Damaging."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            in_silico_predictions=InSilicoPredictions(
                sift_score=0.0,
                polyphen2_score=1.0,
                revel_score=0.95,
                cadd_phred=33.0,
                consensus="Damaging",
                damaging_count=4,
                benign_count=0,
                total_predictors=4,
            ),
        )
        result = self.scorer.score_variant("BRAF", "V600E", evidence_bundle=bundle)
        codes = [c.code for c in result.applied_codes]
        assert "OP1" in codes
        op1 = next(c for c in result.applied_codes if c.code == "OP1")
        assert "4/4" in op1.evidence

    def test_auto_detect_sbp1_with_benign_predictions(self):
        """SBP1 should fire when in-silico consensus is Benign."""
        bundle = EvidenceBundle(
            gene="SOMEGENE",
            variant="A100T",
            in_silico_predictions=InSilicoPredictions(
                sift_score=0.5,
                polyphen2_score=0.05,
                revel_score=0.1,
                consensus="Benign",
                damaging_count=0,
                benign_count=3,
                total_predictors=3,
            ),
        )
        result = self.scorer.score_variant("SOMEGENE", "A100T", evidence_bundle=bundle)
        codes = [c.code for c in result.applied_codes]
        assert "SBP1" in codes

    def test_om1_fallback_without_domain_data(self):
        """OM1 should fall back to gene-level classification when no domain data."""
        bundle = EvidenceBundle(gene="BRAF", variant="V600E")
        result = self.scorer.score_variant("BRAF", "V600E", evidence_bundle=bundle)
        codes = [c.code for c in result.applied_codes]
        assert "OM1" in codes
        om1 = next(c for c in result.applied_codes if c.code == "OM1")
        assert "gene-level" in om1.evidence


class TestACMGAMPHelper:
    def setup_method(self):
        self.helper = ACMGAMPHelper()

    def test_explain_pathogenic(self):
        """Should explain a Pathogenic classification."""
        cv = ClinVarVariant(
            variation_id=1,
            clinical_significance="Pathogenic",
            review_status="reviewed by expert panel",
            review_stars="★★★",
            submitter_count=10,
            conflicting=False,
            submitter_classifications=[
                {"submitter": "Lab1", "classification": "Pathogenic"},
                {"submitter": "Lab2", "classification": "Pathogenic"},
            ],
        )
        result = self.helper.explain_clinvar_classification(cv)
        assert result.aggregate_classification == "Pathogenic"
        assert result.review_stars == "★★★"
        assert not result.conflicting_interpretations
        assert "disease-causing" in result.explanation.lower()

    def test_explain_conflicting(self):
        """Should flag conflicting interpretations."""
        cv = ClinVarVariant(
            variation_id=1,
            clinical_significance="Conflicting interpretations",
            review_status="criteria provided, single submitter",
            submitter_count=3,
            conflicting=True,
            submitter_classifications=[
                {"submitter": "Lab1", "classification": "Pathogenic"},
                {"submitter": "Lab2", "classification": "Likely benign"},
            ],
        )
        result = self.helper.explain_clinvar_classification(cv)
        assert result.conflicting_interpretations
        assert "CONFLICTING" in result.explanation

    def test_get_criteria_reference_pvs1(self):
        """Should return PVS1 description."""
        ref = self.helper.get_criteria_reference("PVS1")
        assert "PVS1" in ref
        assert "null variant" in ref.lower()

    def test_get_criteria_reference_unknown(self):
        """Unknown code should return helpful error."""
        ref = self.helper.get_criteria_reference("FAKE1")
        assert "Unknown" in ref

    def test_get_all_criteria(self):
        """Should return a complete reference."""
        all_ref = self.helper.get_all_criteria()
        assert "PVS1" in all_ref
        assert "BA1" in all_ref
        assert "BP7" in all_ref
        assert "Combining Criteria" in all_ref
