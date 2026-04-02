"""Tests for report and table formatters."""

from variant_mcp.formatters.reports import ReportFormatter
from variant_mcp.formatters.tables import TableFormatter
from variant_mcp.models.classification import (
    AMPTierResult,
    AppliedEvidenceCode,
    OncogenicityScoringResult,
    VariantClassificationReport,
)
from variant_mcp.models.evidence import (
    CIViCEvidenceItem,
    ClinVarVariant,
    EvidenceBundle,
    InSilicoPredictions,
    OncoKBAnnotation,
    ProteinDomain,
)


class TestReportFormatter:
    def test_format_evidence_report_empty(self):
        """Empty bundle should still produce valid report."""
        bundle = EvidenceBundle(gene="BRAF", variant="V600E")
        report = ReportFormatter.format_evidence_report(bundle)
        assert "BRAF" in report
        assert "V600E" in report
        assert "DISCLAIMER" in report or "RESEARCH" in report

    def test_format_evidence_report_with_civic(self):
        """Report with CIViC data should include evidence items."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            civic_evidence=[
                CIViCEvidenceItem(
                    id=100,
                    evidence_type="PREDICTIVE",
                    evidence_level="A",
                    significance="SENSITIVITY",
                    disease="Melanoma",
                    therapies=["Vemurafenib"],
                ),
            ],
        )
        report = ReportFormatter.format_evidence_report(bundle)
        assert "EID100" in report
        assert "Melanoma" in report
        assert "Vemurafenib" in report

    def test_format_classification_report_somatic(self):
        """Somatic classification report should include AMP tier."""
        report_data = VariantClassificationReport(
            gene="BRAF",
            variant="V600E",
            variant_origin="somatic",
            amp_tier=AMPTierResult(
                tier="Tier I",
                tier_name="Tier I — Strong Clinical Significance",
                evidence_level="Level A",
                confidence="HIGH",
                evidence_trail=["OncoKB Level 1"],
            ),
            oncogenicity=OncogenicityScoringResult(
                total_points=14,
                classification="Oncogenic",
                applied_codes=[
                    AppliedEvidenceCode(code="OS1", points=4, description="Test"),
                ],
            ),
            sources_queried=["CIViC", "ClinVar", "OncoKB"],
        )
        report = ReportFormatter.format_classification_report(report_data)
        assert "Tier I" in report
        assert "Oncogenic" in report
        assert "AMP/ASCO/CAP" in report
        assert "ClinGen" in report

    def test_format_amp_tier(self):
        """AMP tier formatting."""
        result = AMPTierResult(
            tier="Tier II",
            tier_name="Tier II — Potential Clinical Significance",
            evidence_level="Level C",
            confidence="MODERATE",
            evidence_trail=["Some evidence"],
            sources_used=["CIViC"],
        )
        text = ReportFormatter.format_amp_tier(result)
        assert "Tier II" in text
        assert "Level C" in text

    def test_format_pathogenicity_summary(self):
        """Pathogenicity summary should include both sections."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            oncokb_annotation=OncoKBAnnotation(
                oncogenic="Oncogenic",
                known_effect="Gain-of-function",
            ),
        )
        report = ReportFormatter.format_pathogenicity_summary(bundle)
        assert "Pathogenic" in report or "Oncogenic" in report
        assert "Benign" in report

    def test_format_source_comparison(self):
        """Source comparison should produce a table."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            civic_evidence=[
                CIViCEvidenceItem(id=1, evidence_level="A", significance="SENSITIVITY"),
            ],
            clinvar_variants=[
                ClinVarVariant(
                    variation_id=1,
                    clinical_significance="Pathogenic",
                    review_stars="★★★",
                ),
            ],
        )
        report = ReportFormatter.format_source_comparison(bundle)
        assert "CIViC" in report
        assert "ClinVar" in report
        assert "OncoKB" in report


    def test_format_evidence_report_with_domains(self):
        """Evidence report should show UniProt protein domains."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            protein_domains=[
                ProteinDomain(
                    name="Protein kinase",
                    domain_type="Domain",
                    start_pos=457,
                    end_pos=717,
                ),
            ],
        )
        report = ReportFormatter.format_evidence_report(bundle)
        assert "UniProt Protein Domains" in report
        assert "Protein kinase" in report
        assert "457" in report
        assert "717" in report

    def test_format_evidence_report_with_predictions(self):
        """Evidence report should show in-silico prediction table."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            in_silico_predictions=InSilicoPredictions(
                sift_score=0.0,
                sift_prediction="deleterious",
                polyphen2_score=1.0,
                polyphen2_prediction="probably_damaging",
                revel_score=0.95,
                cadd_phred=33.0,
                consensus="Damaging",
                damaging_count=4,
                total_predictors=4,
            ),
        )
        report = ReportFormatter.format_evidence_report(bundle)
        assert "In-Silico Predictions" in report
        assert "SIFT" in report
        assert "PolyPhen-2" in report
        assert "REVEL" in report
        assert "CADD" in report
        assert "Damaging" in report
        assert "4/4" in report

    def test_format_pathogenicity_summary_with_domains_and_predictions(self):
        """Pathogenicity summary should include domain and prediction context."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            protein_domains=[
                ProteinDomain(name="Protein kinase", domain_type="Domain", start_pos=457, end_pos=717),
            ],
            in_silico_predictions=InSilicoPredictions(
                consensus="Damaging",
                damaging_count=5,
                total_predictors=5,
            ),
        )
        report = ReportFormatter.format_pathogenicity_summary(bundle)
        assert "UniProt" in report
        assert "Protein kinase" in report
        assert "OM1" in report
        assert "5/5" in report
        assert "OP1" in report or "PP3" in report

    def test_format_pathogenicity_summary_benign_predictions(self):
        """Pathogenicity summary should show benign prediction evidence."""
        bundle = EvidenceBundle(
            gene="SOMEGENE",
            variant="A100T",
            in_silico_predictions=InSilicoPredictions(
                consensus="Benign",
                benign_count=3,
                total_predictors=3,
            ),
        )
        report = ReportFormatter.format_pathogenicity_summary(bundle)
        assert "Benign" in report
        assert "3/3" in report
        assert "SBP1" in report or "BP4" in report

    def test_format_source_comparison_with_domains(self):
        """Source comparison should show protein domain context."""
        bundle = EvidenceBundle(
            gene="BRAF",
            variant="V600E",
            protein_domains=[
                ProteinDomain(name="Protein kinase", domain_type="Domain", start_pos=457, end_pos=717),
            ],
            in_silico_predictions=InSilicoPredictions(
                consensus="Damaging",
                damaging_count=4,
                benign_count=0,
                total_predictors=4,
            ),
        )
        report = ReportFormatter.format_source_comparison(bundle)
        assert "Protein Domain Context" in report
        assert "Protein kinase" in report
        assert "In-Silico Prediction Consensus" in report
        assert "Damaging" in report


class TestTableFormatter:
    def test_frameworks_reference(self):
        """Should include all three frameworks."""
        ref = TableFormatter.format_frameworks_reference()
        assert "AMP/ASCO/CAP" in ref
        assert "ClinGen/CGC/VICC" in ref
        assert "ACMG/AMP" in ref
        assert "When to Use" in ref

    def test_evidence_levels_table(self):
        """Should include all CIViC evidence levels."""
        table = TableFormatter.format_evidence_levels_table()
        assert "Level" in table
        assert "Validated" in table
        assert "Preclinical" in table

    def test_oncokb_levels_table(self):
        """Should include OncoKB levels."""
        table = TableFormatter.format_oncokb_levels_table()
        assert "LEVEL_1" in table
        assert "FDA" in table
