"""Markdown report generators for variant classification and evidence."""

from __future__ import annotations

from datetime import UTC, datetime

from variant_mcp.constants import (
    CIVIC_EVIDENCE_LEVELS,
    DISCLAIMER,
    ONCOKB_LEVELS,
    REFERENCE_PMIDS,
)
from variant_mcp.models.classification import (
    ACMGInterpretation,
    AMPTierResult,
    OncogenicityScoringResult,
    VariantClassificationReport,
)
from variant_mcp.models.evidence import EvidenceBundle


class ReportFormatter:
    """Generates formatted Markdown reports for variant data."""

    @staticmethod
    def format_classification_report(report: VariantClassificationReport) -> str:
        """Format a full variant classification report."""
        lines: list[str] = []
        lines.append("# Variant Classification Report\n")
        lines.append(
            f"**Gene**: {report.gene} | **Variant**: {report.variant} | "
            f"**Origin**: {report.variant_origin.capitalize()}"
        )
        if report.disease:
            lines.append(f"**Disease Context**: {report.disease}")
        lines.append(f"\n**Sources queried**: {', '.join(report.sources_queried)}")
        lines.append(f"**Retrieved**: {datetime.now(UTC).strftime('%Y-%m-%d %H:%M UTC')}")
        lines.append("\n---\n")

        # AMP/ASCO/CAP Classification
        if report.amp_tier:
            lines.append(ReportFormatter.format_amp_tier(report.amp_tier))
            lines.append("\n---\n")

        # Oncogenicity SOP
        if report.oncogenicity:
            lines.append(ReportFormatter.format_oncogenicity(report.oncogenicity))
            lines.append("\n---\n")

        # ACMG Interpretation
        if report.acmg_interpretation:
            lines.append(ReportFormatter.format_acmg_interpretation(report.acmg_interpretation))
            lines.append("\n---\n")

        # Citations
        lines.append("## Citations\n")
        if report.amp_tier:
            lines.append(
                f'1. Li et al. (2017). "Standards and Guidelines for the Interpretation '
                f'and Reporting..." J Mol Diagn. PMID: {REFERENCE_PMIDS["AMP_ASCO_CAP"]}'
            )
        if report.oncogenicity:
            lines.append(
                f'2. Horak et al. (2022). "Standards for the classification of '
                f'pathogenicity of somatic variants..." Genet Med. PMID: {REFERENCE_PMIDS["ONCOGENICITY_SOP"]}'
            )
        if report.acmg_interpretation:
            lines.append(
                f'3. Richards et al. (2015). "Standards and guidelines for the '
                f'interpretation of sequence variants..." Genet Med. PMID: {REFERENCE_PMIDS["ACMG_AMP"]}'
            )

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)

    @staticmethod
    def format_amp_tier(result: AMPTierResult) -> str:
        """Format AMP/ASCO/CAP tier result."""
        lines: list[str] = []
        lines.append("## AMP/ASCO/CAP Classification\n")
        level_str = f" ({result.evidence_level})" if result.evidence_level else ""
        lines.append(f"**{result.tier_name}{level_str}**")
        lines.append(f"**Confidence**: {result.confidence}\n")

        if result.evidence_trail:
            lines.append("**Evidence trail**:")
            for item in result.evidence_trail:
                lines.append(f"- {item}")

        if result.sources_used:
            lines.append(f"\n*Sources*: {', '.join(result.sources_used)}")
        lines.append(f"*Reference*: {result.guideline_reference}")
        return "\n".join(lines)

    @staticmethod
    def format_oncogenicity(result: OncogenicityScoringResult) -> str:
        """Format oncogenicity scoring result."""
        lines: list[str] = []
        lines.append("## ClinGen/CGC/VICC Oncogenicity Assessment\n")
        lines.append(
            f"**Classification: {result.classification.upper()} "
            f"(Score: {result.total_points} points)**"
        )
        if result.auto_detected:
            lines.append("*Evidence codes were auto-detected from available data*")
        lines.append(f"**Confidence**: {result.confidence}\n")

        if result.applied_codes:
            lines.append("| Code | Points | Evidence |")
            lines.append("|------|--------|----------|")
            for code in result.applied_codes:
                lines.append(f"| {code.code} | {code.points:+d} | {code.description} |")
            lines.append(
                f"| **Total** | **{result.total_points:+d}** | "
                f"**{result.classification} (threshold: >=10 Oncogenic, "
                f"6-9 Likely Oncogenic)** |"
            )

        lines.append(f"\n*Reference*: {result.guideline_reference}")
        return "\n".join(lines)

    @staticmethod
    def format_acmg_interpretation(result: ACMGInterpretation) -> str:
        """Format ACMG/AMP interpretation."""
        lines: list[str] = []
        lines.append("## ACMG/AMP Germline Interpretation\n")
        lines.append(result.explanation)
        lines.append(f"\n*Reference*: {result.guideline_reference}")
        return "\n".join(lines)

    @staticmethod
    def format_evidence_report(bundle: EvidenceBundle) -> str:
        """Format a multi-source evidence search report."""
        lines: list[str] = []
        lines.append("# Multi-Source Evidence Report\n")
        lines.append(f"**Gene**: {bundle.gene} | **Variant**: {bundle.variant or 'Any'}")
        if bundle.disease:
            lines.append(f"**Disease**: {bundle.disease}")
        lines.append(f"**Sources with data**: {', '.join(bundle.sources_queried) or 'None'}")
        lines.append(f"**Retrieved**: {datetime.now(UTC).strftime('%Y-%m-%d %H:%M UTC')}")
        lines.append("\n---\n")

        # CIViC section
        if bundle.civic_evidence:
            lines.append("## CIViC (civicdb.org)\n")
            lines.append(f"**Total evidence items**: {len(bundle.civic_evidence)}")

            type_counts: dict[str, int] = {}
            level_counts: dict[str, int] = {}
            for item in bundle.civic_evidence:
                et = item.evidence_type or "Unknown"
                el = item.evidence_level or "?"
                type_counts[et] = type_counts.get(et, 0) + 1
                level_counts[el] = level_counts.get(el, 0) + 1

            lines.append(
                "**Evidence types**: "
                + ", ".join(f"{t} ({c})" for t, c in sorted(type_counts.items()))
            )
            lines.append(
                "**Evidence levels**: "
                + ", ".join(f"Level {lvl} ({cnt})" for lvl, cnt in sorted(level_counts.items()))
            )

            # Top evidence items
            lines.append("\n**Top evidence items**:")
            for item in bundle.civic_evidence[:10]:
                level_desc = CIVIC_EVIDENCE_LEVELS.get(item.evidence_level or "", "")
                therapies_str = ", ".join(item.therapies) if item.therapies else "N/A"
                lines.append(
                    f"- **EID{item.id}** [{item.evidence_type}] "
                    f"Level {item.evidence_level} ({level_desc}) — "
                    f"{item.significance or 'N/A'} | "
                    f"Disease: {item.disease or 'N/A'} | "
                    f"Therapies: {therapies_str}"
                )
                if item.source and item.source.pmid:
                    lines.append(f"  PMID: {item.source.pmid}")
            lines.append("")

        # CIViC Assertions
        if bundle.civic_assertions:
            lines.append("### CIViC Assertions\n")
            for assertion in bundle.civic_assertions:
                lines.append(
                    f"- **AID{assertion.id}**: {assertion.assertion_type} — "
                    f"{assertion.significance} | "
                    f"AMP Level: {assertion.amp_level or 'N/A'} | "
                    f"Disease: {assertion.disease or 'N/A'}"
                )
            lines.append("")

        # ClinVar section
        if bundle.clinvar_variants:
            lines.append("## ClinVar (ncbi.nlm.nih.gov/clinvar)\n")
            for cv in bundle.clinvar_variants:
                lines.append(f"**Variation ID**: {cv.variation_id}")
                lines.append(
                    f"**Classification**: {cv.clinical_significance} "
                    f"({cv.review_stars} — {cv.review_status})"
                )
                if cv.submitter_count:
                    lines.append(f"**Submitters**: {cv.submitter_count}")
                if cv.conflicting:
                    lines.append("**CONFLICTING** interpretations among submitters")
                if cv.conditions:
                    lines.append(f"**Conditions**: {', '.join(cv.conditions[:5])}")
                if cv.clinvar_url:
                    lines.append(f"**URL**: {cv.clinvar_url}")
                lines.append("")

        # OncoKB section
        if bundle.has_oncokb_data:
            ann = bundle.oncokb_annotation
            assert ann is not None
            lines.append("## OncoKB (oncokb.org)\n")
            lines.append(f"**Oncogenic**: {ann.oncogenic or 'N/A'}")
            lines.append(f"**Mutation Effect**: {ann.known_effect or 'N/A'}")
            if ann.highest_sensitive_level:
                level_desc = ONCOKB_LEVELS.get(ann.highest_sensitive_level, "")
                lines.append(
                    f"**Highest Sensitive Level**: {ann.highest_sensitive_level} ({level_desc})"
                )
            if ann.highest_resistance_level:
                level_desc = ONCOKB_LEVELS.get(ann.highest_resistance_level, "")
                lines.append(
                    f"**Highest Resistance Level**: {ann.highest_resistance_level} ({level_desc})"
                )
            if ann.treatments:
                lines.append("**Treatments**:")
                for tx in ann.treatments[:5]:
                    drugs = ", ".join(tx.get("drugs", []))
                    lines.append(f"  - {tx.get('level', '')}: {drugs}")
            if ann.oncokb_url:
                lines.append(f"**URL**: {ann.oncokb_url}")
            lines.append("")

        # MetaKB section
        if bundle.metakb_interpretations:
            lines.append("## VICC MetaKB (search.cancervariants.org)\n")
            for interp in bundle.metakb_interpretations[:10]:
                drugs_str = ", ".join(interp.drugs) if interp.drugs else "N/A"
                lines.append(
                    f"- **{interp.source or 'Unknown'}**: "
                    f"{interp.clinical_significance or 'N/A'} | "
                    f"Disease: {interp.disease or 'N/A'} | "
                    f"Drugs: {drugs_str}"
                )
            lines.append("")

        # Errors
        if bundle.errors:
            lines.append("## Data Source Errors\n")
            for source, error in bundle.errors.items():
                lines.append(f"- **{source}**: {error}")
            lines.append("")

        lines.append(f"---\n\n{DISCLAIMER}")
        return "\n".join(lines)

    @staticmethod
    def format_pathogenicity_summary(
        bundle: EvidenceBundle,
        oncogenicity: OncogenicityScoringResult | None = None,
    ) -> str:
        """Format a pathogenic vs benign evidence summary."""
        lines: list[str] = []
        lines.append(f"# Pathogenicity Summary: {bundle.gene} {bundle.variant}\n")

        # Pathogenic / Oncogenic Evidence
        lines.append("## Pathogenic/Oncogenic Evidence\n")
        has_pathogenic = False

        if bundle.has_clinvar_data:
            for cv in bundle.clinvar_variants:
                sig = (cv.clinical_significance or "").lower()
                if "pathogenic" in sig or "oncogenic" in sig:
                    breakdown_parts = []
                    for sub in cv.submitter_classifications:
                        cls = sub.get("classification", "")
                        if "pathogenic" in cls.lower():
                            breakdown_parts.append(cls)
                    pathogenic_count = len(breakdown_parts)
                    total = cv.submitter_count
                    lines.append(
                        f"- **ClinVar**: {cv.clinical_significance} "
                        f"({cv.review_stars}) — {pathogenic_count}/{total} "
                        f"submitters support pathogenic classification"
                    )
                    has_pathogenic = True

        if bundle.civic_evidence:
            oncogenic_items = [
                e
                for e in bundle.civic_evidence
                if (e.significance or "").lower()
                in ("sensitivity", "oncogenic", "gain of function", "pathogenic")
            ]
            if oncogenic_items:
                level_breakdown: dict[str, int] = {}
                for item in oncogenic_items:
                    lvl = item.evidence_level or "?"
                    level_breakdown[lvl] = level_breakdown.get(lvl, 0) + 1
                breakdown_str = ", ".join(
                    f"{cnt} Level {lvl}" for lvl, cnt in sorted(level_breakdown.items())
                )
                lines.append(
                    f"- **CIViC**: {len(oncogenic_items)} evidence items support "
                    f"oncogenic/pathogenic function ({breakdown_str})"
                )
                has_pathogenic = True

        if bundle.has_oncokb_data:
            ann = bundle.oncokb_annotation
            assert ann is not None
            oncogenic = (ann.oncogenic or "").lower()
            if "oncogenic" in oncogenic:
                lines.append(f"- **OncoKB**: {ann.oncogenic}, {ann.known_effect or 'N/A'}")
                has_pathogenic = True

        if oncogenicity and oncogenicity.total_points >= 6:
            lines.append(
                f"- **Oncogenicity SOP**: {oncogenicity.total_points} points "
                f"-> {oncogenicity.classification}"
            )
            has_pathogenic = True

        if not has_pathogenic:
            lines.append("- No pathogenic/oncogenic evidence found")

        # Benign Evidence
        lines.append("\n## Benign Evidence\n")
        has_benign = False

        if bundle.has_clinvar_data:
            for cv in bundle.clinvar_variants:
                sig = (cv.clinical_significance or "").lower()
                if "benign" in sig:
                    lines.append(f"- **ClinVar**: {cv.clinical_significance} ({cv.review_stars})")
                    has_benign = True

        if bundle.civic_evidence:
            benign_items = [
                e
                for e in bundle.civic_evidence
                if (e.significance or "").lower() in ("benign", "likely benign", "neutral")
            ]
            if benign_items:
                lines.append(
                    f"- **CIViC**: {len(benign_items)} evidence items with benign significance"
                )
                has_benign = True

        if bundle.has_oncokb_data:
            ann = bundle.oncokb_annotation
            assert ann is not None
            if (ann.oncogenic or "").lower() in ("likely neutral", "inconclusive"):
                lines.append(f"- **OncoKB**: {ann.oncogenic}")
                has_benign = True

        if not has_benign:
            lines.append("- No benign evidence found")

        # Consensus Assessment
        lines.append("\n## Consensus Assessment\n")
        if has_pathogenic and not has_benign:
            sources = len(bundle.sources_queried)
            if sources >= 3:
                lines.append(
                    "**STRONG PATHOGENIC/ONCOGENIC** — Multiple independent sources "
                    "agree. Confidence: HIGH"
                )
            elif sources >= 2:
                lines.append(
                    "**LIKELY PATHOGENIC/ONCOGENIC** — Multiple sources support. "
                    "Confidence: MODERATE"
                )
            else:
                lines.append(
                    "**PATHOGENIC/ONCOGENIC EVIDENCE** — Limited sources queried. Confidence: LOW"
                )
        elif has_benign and not has_pathogenic:
            lines.append("**BENIGN** — Evidence supports benign classification.")
        elif has_pathogenic and has_benign:
            lines.append(
                "**CONFLICTING** — Both pathogenic and benign evidence exist. "
                "Manual expert review recommended."
            )
        else:
            lines.append(
                "**INSUFFICIENT DATA** — No strong evidence in either direction. Classified as VUS."
            )

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)

    @staticmethod
    def format_source_comparison(bundle: EvidenceBundle) -> str:
        """Format a cross-database comparison for a variant."""
        lines: list[str] = []
        lines.append(f"# Cross-Database Comparison: {bundle.gene} {bundle.variant}\n")
        lines.append(f"**Retrieved**: {datetime.now(UTC).strftime('%Y-%m-%d %H:%M UTC')}\n")

        lines.append("| Attribute | CIViC | ClinVar | OncoKB |")
        lines.append("|-----------|-------|---------|--------|")

        # Data availability
        civic_status = (
            f"{len(bundle.civic_evidence)} evidence items" if bundle.civic_evidence else "No data"
        )
        clinvar_status = (
            bundle.clinvar_variants[0].clinical_significance
            if bundle.clinvar_variants
            else "No data"
        )
        oncokb_status: str | None = (
            bundle.oncokb_annotation.oncogenic
            if bundle.has_oncokb_data and bundle.oncokb_annotation
            else "No data / No token"
        )
        lines.append(f"| Status | {civic_status} | {clinvar_status} | {oncokb_status} |")

        # Classification
        civic_sig = set()
        for e in bundle.civic_evidence:
            if e.significance:
                civic_sig.add(e.significance)
        civic_class = ", ".join(civic_sig) if civic_sig else "N/A"

        clinvar_class = "N/A"
        if bundle.clinvar_variants:
            clinvar_class = bundle.clinvar_variants[0].clinical_significance or "N/A"

        oncokb_class = "N/A"
        if bundle.has_oncokb_data and bundle.oncokb_annotation:
            oncokb_class = bundle.oncokb_annotation.oncogenic or "N/A"

        lines.append(f"| Classification | {civic_class} | {clinvar_class} | {oncokb_class} |")

        # Evidence quality
        civic_quality = "N/A"
        if bundle.civic_evidence:
            levels = [e.evidence_level for e in bundle.civic_evidence if e.evidence_level]
            if levels:
                best = min(levels)
                civic_quality = f"Best: Level {best}"

        clinvar_quality = "N/A"
        if bundle.clinvar_variants:
            clinvar_quality = bundle.clinvar_variants[0].review_stars or "N/A"

        oncokb_quality = "N/A"
        if (
            bundle.has_oncokb_data
            and bundle.oncokb_annotation
            and bundle.oncokb_annotation.highest_sensitive_level
        ):
            oncokb_quality = bundle.oncokb_annotation.highest_sensitive_level

        lines.append(
            f"| Evidence Quality | {civic_quality} | {clinvar_quality} | {oncokb_quality} |"
        )

        # Therapeutic info
        civic_therapies: set[str] = set()
        for e in bundle.civic_evidence:
            for t in e.therapies:
                civic_therapies.add(t)
        civic_tx = ", ".join(list(civic_therapies)[:3]) if civic_therapies else "N/A"

        oncokb_tx = "N/A"
        if (
            bundle.has_oncokb_data
            and bundle.oncokb_annotation
            and bundle.oncokb_annotation.treatments
        ):
            drugs = []
            for tx in bundle.oncokb_annotation.treatments[:3]:
                drugs.extend(tx.get("drugs", []))
            oncokb_tx = ", ".join(drugs[:3]) if drugs else "N/A"

        lines.append(f"| Therapies | {civic_tx} | N/A | {oncokb_tx} |")

        # Concordance assessment
        lines.append("\n## Concordance Assessment\n")
        classifications = []
        if civic_sig:
            classifications.append(("CIViC", civic_class))
        if clinvar_class != "N/A":
            classifications.append(("ClinVar", clinvar_class))
        if oncokb_class != "N/A":
            classifications.append(("OncoKB", oncokb_class))

        if len(classifications) >= 2:
            lines.append(
                "**Sources compared**: " + ", ".join(f"{s}: {c}" for s, c in classifications)
            )
        else:
            lines.append("Insufficient data from multiple sources for comparison.")

        if bundle.errors:
            lines.append("\n## Errors\n")
            for source, error in bundle.errors.items():
                lines.append(f"- **{source}**: {error}")

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)
