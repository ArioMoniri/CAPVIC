"""AMP/ASCO/CAP 4-Tier Somatic Variant Classifier.

Reference: Li et al., J Mol Diagn 2017. PMID: 27993330
"""

from __future__ import annotations

import logging

from variant_mcp.constants import (
    AMP_TIER_DEFINITIONS,
    ONCOKB_LEVELS,
)
from variant_mcp.models.classification import AMPTierResult
from variant_mcp.models.evidence import EvidenceBundle

logger = logging.getLogger(__name__)


class AMPTierClassifier:
    """AMP/ASCO/CAP somatic variant clinical significance tiering.

    Tier I: Strong clinical significance
      Level A: FDA-approved therapy, guideline-recommended
      Level B: Well-powered clinical studies with consensus
    Tier II: Potential clinical significance
      Level C: FDA-approved for different tumor type, or investigational in trials
      Level D: Preclinical/case report evidence
    Tier III: Unknown clinical significance
    Tier IV: Benign or Likely Benign
    """

    def classify(self, evidence_bundle: EvidenceBundle) -> AMPTierResult:
        """Classify a somatic variant using aggregated evidence.

        Logic:
        1. Check OncoKB levels: L1 -> Tier I/A, L2 -> Tier I/B, L3A -> Tier II/C, etc.
        2. Check CIViC evidence levels and assertions
        3. Cross-reference ClinVar for benign signals
        4. Population frequency checks (Tier IV support)
        5. No evidence -> Tier III
        """
        trail: list[str] = []
        sources: list[str] = []
        tier = "III"
        level: str | None = None
        confidence = "LOW"

        # 1. Check OncoKB
        oncokb_tier, oncokb_level, oncokb_trail = self._check_oncokb(evidence_bundle)
        if oncokb_tier:
            tier = oncokb_tier
            level = oncokb_level
            trail.extend(oncokb_trail)
            sources.append("OncoKB")

        # 2. Check CIViC assertions (higher confidence than evidence items)
        civic_a_tier, civic_a_level, civic_a_trail = self._check_civic_assertions(evidence_bundle)
        if civic_a_tier and self._tier_rank(civic_a_tier) < self._tier_rank(tier):
            tier = civic_a_tier
            level = civic_a_level
            trail.extend(civic_a_trail)
        elif civic_a_trail:
            trail.extend(civic_a_trail)
        if evidence_bundle.civic_assertions:
            sources.append("CIViC Assertions")

        # 3. Check CIViC evidence items
        civic_e_tier, civic_e_level, civic_e_trail = self._check_civic_evidence(evidence_bundle)
        if civic_e_tier and self._tier_rank(civic_e_tier) < self._tier_rank(tier):
            tier = civic_e_tier
            level = civic_e_level
            trail.extend(civic_e_trail)
        elif civic_e_trail:
            trail.extend(civic_e_trail)
        if evidence_bundle.civic_evidence:
            sources.append("CIViC Evidence")

        # 4. Check ClinVar for benign/pathogenic signals
        clinvar_tier, clinvar_trail = self._check_clinvar(evidence_bundle)
        if clinvar_tier == "IV" and tier == "III":
            tier = "IV"
            trail.extend(clinvar_trail)
        elif clinvar_trail:
            trail.extend(clinvar_trail)
        if evidence_bundle.has_clinvar_data:
            sources.append("ClinVar")

        # Determine confidence
        confidence = self._assess_confidence(tier, level, trail, sources)

        tier_name = AMP_TIER_DEFINITIONS.get(tier, {}).get("name", f"Tier {tier}")

        return AMPTierResult(
            tier=f"Tier {tier}",
            tier_name=tier_name,
            evidence_level=f"Level {level}" if level else None,
            confidence=confidence,
            evidence_trail=trail,
            sources_used=sources,
        )

    def _check_oncokb(self, bundle: EvidenceBundle) -> tuple[str | None, str | None, list[str]]:
        """Check OncoKB therapeutic levels for AMP tier mapping."""
        if not bundle.has_oncokb_data:
            return None, None, []

        annotation = bundle.oncokb_annotation
        if annotation is None:
            return None, None, []
        trail: list[str] = []
        tier: str | None = None
        level: str | None = None

        sensitive = annotation.highest_sensitive_level or ""
        resistance = annotation.highest_resistance_level or ""

        # Map OncoKB levels to AMP tiers
        if sensitive in ("LEVEL_1",):
            tier, level = "I", "A"
            trail.append(f"OncoKB Level 1: {ONCOKB_LEVELS.get(sensitive, sensitive)}")
        elif sensitive in ("LEVEL_2",):
            tier, level = "I", "B"
            trail.append(f"OncoKB Level 2: {ONCOKB_LEVELS.get(sensitive, sensitive)}")
        elif sensitive in ("LEVEL_3A",):
            tier, level = "II", "C"
            trail.append(f"OncoKB Level 3A: {ONCOKB_LEVELS.get(sensitive, sensitive)}")
        elif sensitive in ("LEVEL_3B", "LEVEL_4"):
            tier, level = "II", "D"
            trail.append(f"OncoKB {sensitive}: {ONCOKB_LEVELS.get(sensitive, sensitive)}")

        if resistance in ("LEVEL_R1",):
            trail.append(f"OncoKB R1 resistance: {ONCOKB_LEVELS.get(resistance, resistance)}")
        elif resistance in ("LEVEL_R2",):
            trail.append(f"OncoKB R2 resistance: {ONCOKB_LEVELS.get(resistance, resistance)}")

        # Oncogenicity status
        if annotation.oncogenic:
            trail.append(f"OncoKB Oncogenicity: {annotation.oncogenic}")
        if annotation.known_effect:
            trail.append(f"OncoKB Mutation Effect: {annotation.known_effect}")

        # Treatments
        for tx in annotation.treatments[:3]:
            drugs = ", ".join(tx.get("drugs", []))
            tx_level = tx.get("level", "")
            if drugs:
                trail.append(f"OncoKB Treatment ({tx_level}): {drugs}")

        return tier, level, trail

    def _check_civic_assertions(
        self, bundle: EvidenceBundle
    ) -> tuple[str | None, str | None, list[str]]:
        """Check CIViC assertions for AMP tier signals."""
        if not bundle.civic_assertions:
            return None, None, []

        trail: list[str] = []
        best_tier: str | None = None
        best_level: str | None = None

        for assertion in bundle.civic_assertions:
            amp = assertion.amp_level
            if amp:
                amp_upper = amp.upper().replace(" ", "")
                a_tier, a_level = self._parse_amp_level(amp_upper)
                if a_tier and (
                    best_tier is None or self._tier_rank(a_tier) < self._tier_rank(best_tier)
                ):
                    best_tier = a_tier
                    best_level = a_level

                trail.append(
                    f"CIViC Assertion {assertion.id}: {assertion.assertion_type} — "
                    f"AMP Level {amp}, significance={assertion.significance}"
                )
            if assertion.nccn_guideline:
                trail.append(
                    f"CIViC Assertion {assertion.id}: NCCN guideline — {assertion.nccn_guideline}"
                )

        return best_tier, best_level, trail

    def _check_civic_evidence(
        self, bundle: EvidenceBundle
    ) -> tuple[str | None, str | None, list[str]]:
        """Check CIViC evidence items for AMP tier signals."""
        if not bundle.civic_evidence:
            return None, None, []

        trail: list[str] = []
        best_tier: str | None = None
        best_level: str | None = None

        # Count by evidence level
        level_counts: dict[str, int] = {}
        for item in bundle.civic_evidence:
            el = item.evidence_level or ""
            level_counts[el] = level_counts.get(el, 0) + 1

        # Map CIViC evidence levels to AMP tiers
        # Level A (Validated) -> Tier I/A
        # Level B (Clinical) -> Tier I/B or Tier II/C
        # Level C (Case Study) -> Tier II/D
        # Level D (Preclinical) -> Tier II/D
        if level_counts.get("A", 0) > 0:
            best_tier, best_level = "I", "A"
            trail.append(f"CIViC: {level_counts['A']} Level A (Validated) evidence items")
        elif level_counts.get("B", 0) > 0:
            best_tier, best_level = "I", "B"
            trail.append(f"CIViC: {level_counts['B']} Level B (Clinical) evidence items")
        elif level_counts.get("C", 0) > 0:
            best_tier, best_level = "II", "D"
            trail.append(f"CIViC: {level_counts['C']} Level C (Case Study) evidence items")
        elif level_counts.get("D", 0) > 0:
            best_tier, best_level = "II", "D"
            trail.append(f"CIViC: {level_counts['D']} Level D (Preclinical) evidence items")

        # Evidence type breakdown
        type_counts: dict[str, int] = {}
        for item in bundle.civic_evidence:
            et = item.evidence_type or "Unknown"
            type_counts[et] = type_counts.get(et, 0) + 1
        for etype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
            trail.append(f"CIViC evidence type: {etype} ({count} items)")

        return best_tier, best_level, trail

    def _check_clinvar(self, bundle: EvidenceBundle) -> tuple[str | None, list[str]]:
        """Check ClinVar for benign/pathogenic signals."""
        if not bundle.has_clinvar_data:
            return None, []

        trail: list[str] = []
        tier: str | None = None

        for cv in bundle.clinvar_variants:
            sig = (cv.clinical_significance or "").lower()
            stars = cv.review_stars or ""
            status = cv.review_status or ""

            trail.append(
                f"ClinVar {cv.variation_id}: {cv.clinical_significance} ({stars} — {status})"
            )

            # Only assign Tier IV if the classification is clearly benign,
            # not when pathogenic/benign conflict exists
            if "pathogenic" not in sig:
                if "benign" in sig and "likely" not in sig:
                    tier = "IV"
                    trail.append("ClinVar classification supports Tier IV (Benign)")
                elif "likely benign" in sig:
                    tier = "IV"
                    trail.append("ClinVar classification supports Tier IV (Likely Benign)")

            if cv.conflicting:
                trail.append("ClinVar: CONFLICTING interpretations among submitters")

        return tier, trail

    @staticmethod
    def _parse_amp_level(amp_str: str) -> tuple[str | None, str | None]:
        """Parse AMP level string like 'TIERI-LEVELA' into (tier, level)."""
        amp_str = amp_str.upper().replace("_", "").replace("-", "").replace(" ", "")
        # Check longest tier strings first to avoid prefix collision
        # (e.g., "TIERI" matching before "TIERII")
        for tier_val in ("IV", "III", "II", "I"):
            for level_val in ("A", "B", "C", "D"):
                if f"TIER{tier_val}" in amp_str and f"LEVEL{level_val}" in amp_str:
                    return tier_val, level_val
            if f"TIER{tier_val}" in amp_str:
                return tier_val, None
        return None, None

    @staticmethod
    def _tier_rank(tier: str) -> int:
        """Lower rank = stronger evidence."""
        return {"I": 1, "II": 2, "III": 3, "IV": 4}.get(tier, 5)

    @staticmethod
    def _assess_confidence(
        tier: str,
        level: str | None,
        trail: list[str],
        sources: list[str],
    ) -> str:
        """Assess confidence based on evidence quantity and concordance."""
        if len(sources) >= 3 and tier in ("I", "II"):
            return "HIGH"
        if len(sources) >= 2 and tier in ("I", "II"):
            return "MODERATE"
        if len(trail) >= 3:
            return "MODERATE"
        return "LOW"
