"""ACMG/AMP Germline Pathogenicity Helper.

Reference: Richards et al., Genet Med 2015. PMID: 25741868
"""

from __future__ import annotations

import logging

from variant_mcp.constants import ACMG_CRITERIA, CLINVAR_REVIEW_STARS
from variant_mcp.models.classification import ACMGInterpretation
from variant_mcp.models.evidence import ClinVarVariant

logger = logging.getLogger(__name__)


class ACMGAMPHelper:
    """Reference helper for ACMG/AMP germline variant pathogenicity guidelines.

    NOTE: This does NOT perform full ACMG classification (which requires
    manual expert review). It helps teams understand ClinVar classifications
    and the evidence criteria framework.
    """

    def explain_clinvar_classification(self, clinvar_data: ClinVarVariant) -> ACMGInterpretation:
        """Explain a ClinVar variant's classification in ACMG/AMP context.

        Covers: aggregate classification, review status stars, conflicting
        interpretations, submitter details, and evaluation date.
        """
        # Build submitter breakdown
        breakdown: dict[str, int] = {}
        for sub in clinvar_data.submitter_classifications:
            cls = sub.get("classification", "Unknown")
            breakdown[cls] = breakdown.get(cls, 0) + 1

        explanation_parts: list[str] = []

        # Classification explanation
        sig = clinvar_data.clinical_significance or "Not provided"
        explanation_parts.append(f"**Aggregate Classification**: {sig}")
        explanation_parts.append(self._explain_significance(sig))

        # Review status
        review = clinvar_data.review_status or "Not provided"
        stars = clinvar_data.review_stars or CLINVAR_REVIEW_STARS.get(review.lower(), "")
        explanation_parts.append(f"\n**Review Status**: {review} {stars}")
        explanation_parts.append(self._explain_review_status(review))

        # Submitter details
        count = clinvar_data.submitter_count
        explanation_parts.append(f"\n**Submitter Count**: {count}")
        if breakdown:
            explanation_parts.append("**Submitter Breakdown**:")
            for cls, num in sorted(breakdown.items(), key=lambda x: -x[1]):
                explanation_parts.append(f"  - {cls}: {num} submitter(s)")

        # Conflicting interpretations
        if clinvar_data.conflicting:
            explanation_parts.append(
                "\n**CONFLICTING INTERPRETATIONS**: Submitters disagree on the "
                "classification. This variant may warrant additional review."
            )
        else:
            explanation_parts.append("\n**No conflicting interpretations** among submitters.")

        # Last evaluated
        if clinvar_data.last_evaluated:
            explanation_parts.append(f"\n**Last Evaluated**: {clinvar_data.last_evaluated}")

        return ACMGInterpretation(
            aggregate_classification=sig,
            review_status=review,
            review_stars=stars,
            submitter_count=count,
            conflicting_interpretations=clinvar_data.conflicting,
            submitter_breakdown=breakdown,
            last_evaluated=clinvar_data.last_evaluated,
            explanation="\n".join(explanation_parts),
        )

    def get_criteria_reference(self, criteria_code: str) -> str:
        """Return description and context for a given ACMG criteria code."""
        code = criteria_code.upper().strip()
        description = ACMG_CRITERIA.get(code)

        if not description:
            available = ", ".join(sorted(ACMG_CRITERIA.keys()))
            return f"Unknown ACMG/AMP criteria code: '{code}'. Available codes: {available}"

        strength = self._get_strength(code)
        direction = "Pathogenic" if code.startswith(("P", "PM", "PP")) else "Benign"
        if code.startswith("BA"):
            direction = "Benign (standalone)"

        return (
            f"**{code}** — {direction} ({strength})\n\n"
            f"{description}\n\n"
            f"*Reference: Richards et al., Genet Med 2015. PMID: 25741868*"
        )

    def get_all_criteria(self) -> str:
        """Return a formatted reference of all ACMG/AMP criteria."""
        sections = {
            "Pathogenic — Very Strong": ["PVS1"],
            "Pathogenic — Strong": ["PS1", "PS2", "PS3", "PS4"],
            "Pathogenic — Moderate": ["PM1", "PM2", "PM3", "PM4", "PM5", "PM6"],
            "Pathogenic — Supporting": ["PP1", "PP2", "PP3", "PP4", "PP5"],
            "Benign — Standalone": ["BA1"],
            "Benign — Strong": ["BS1", "BS2", "BS3", "BS4"],
            "Benign — Supporting": ["BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"],
        }

        lines = ["# ACMG/AMP Germline Variant Classification Criteria\n"]
        lines.append("*Richards et al., Genet Med 2015. PMID: 25741868*\n")

        for section, codes in sections.items():
            lines.append(f"\n## {section}\n")
            for code in codes:
                desc = ACMG_CRITERIA.get(code, "")
                lines.append(f"- **{code}**: {desc}")

        lines.append("\n## Combining Criteria for Classification\n")
        lines.append(
            "- **Pathogenic**: PVS1 + >=1 Strong; OR >=2 Strong; OR 1 Strong + >=3 Moderate/Supporting"
        )
        lines.append(
            "- **Likely Pathogenic**: PVS1 + 1 Moderate; OR 1 Strong + 1-2 Moderate; OR 1 Strong + >=2 Supporting; OR >=3 Moderate; OR 2 Moderate + >=2 Supporting"
        )
        lines.append("- **Likely Benign**: 1 Strong + 1 Supporting benign")
        lines.append("- **Benign**: BA1 standalone; OR >=2 Strong benign")
        lines.append("- **VUS**: Criteria not met for any other classification")

        return "\n".join(lines)

    @staticmethod
    def _explain_significance(sig: str) -> str:
        """Explain what a ClinVar classification means."""
        sig_lower = sig.lower()
        explanations = {
            "pathogenic": (
                "This variant is considered disease-causing with high confidence. "
                "Multiple lines of strong evidence support pathogenicity."
            ),
            "likely pathogenic": (
                "There is >90% certainty that this variant is disease-causing. "
                "Strong evidence supports pathogenicity but may not meet all "
                "criteria for definitive classification."
            ),
            "uncertain significance": (
                "There is insufficient evidence to classify this variant as "
                "pathogenic or benign. Additional research or clinical data may "
                "help resolve the classification."
            ),
            "likely benign": (
                "There is >90% certainty that this variant does NOT cause disease. "
                "Evidence suggests the variant is benign but not fully conclusive."
            ),
            "benign": (
                "This variant is considered NOT disease-causing with high confidence. "
                "Strong evidence supports a benign classification."
            ),
        }
        for key, explanation in explanations.items():
            if key in sig_lower:
                return explanation
        if "conflicting" in sig_lower:
            return (
                "Different submitters have assigned different classifications. "
                "Review individual submitter evidence to understand the basis "
                "for each interpretation."
            )
        return "Classification not recognized or not standard ACMG/AMP terminology."

    @staticmethod
    def _explain_review_status(review: str) -> str:
        """Explain what ClinVar review status means."""
        review_lower = review.lower()
        if "practice guideline" in review_lower:
            return (
                "Highest confidence — reviewed and endorsed by a practice guideline "
                "(e.g., ACMG, AMP, NCCN)."
            )
        if "expert panel" in review_lower:
            return "High confidence — reviewed by an expert panel (e.g., ClinGen, ENIGMA, InSiGHT)."
        if "multiple submitters, no conflicts" in review_lower:
            return (
                "Moderate confidence — multiple independent laboratories agree "
                "on the classification."
            )
        if "single submitter" in review_lower:
            return (
                "Lower confidence — only one laboratory has submitted a "
                "classification with criteria."
            )
        if "no assertion criteria" in review_lower:
            return "Lowest confidence — no standard criteria were provided with the classification."
        return "Review status not recognized."

    @staticmethod
    def _get_strength(code: str) -> str:
        """Get the evidence strength of an ACMG criteria code."""
        if code.startswith("PVS"):
            return "Very Strong"
        if code.startswith("PS"):
            return "Strong"
        if code.startswith("PM"):
            return "Moderate"
        if code.startswith("PP"):
            return "Supporting"
        if code.startswith("BA"):
            return "Standalone"
        if code.startswith("BS"):
            return "Strong"
        if code.startswith("BP"):
            return "Supporting"
        return "Unknown"
