"""ClinGen/CGC/VICC Oncogenicity Standard Operating Procedure Scorer.

Reference: Horak et al., Genet Med 2022. PMID: 35101336
"""

from __future__ import annotations

import logging
import re

from variant_mcp.constants import (
    KNOWN_ONCOGENES,
    KNOWN_TUMOR_SUPPRESSORS,
    NULL_VARIANT_TYPES,
    ONCOGENICITY_EVIDENCE_CODES,
)
from variant_mcp.models.classification import AppliedEvidenceCode, OncogenicityScoringResult
from variant_mcp.models.evidence import EvidenceBundle

logger = logging.getLogger(__name__)


class OncogenicityScorer:
    """ClinGen/CGC/VICC SOP for classification of oncogenicity of somatic variants.

    Uses a point-based system:
        >= 10 points -> Oncogenic
        6-9 points   -> Likely Oncogenic
        0-5 points   -> VUS
        -1 to -5     -> VUS
        -6 to -9     -> Likely Benign
        <= -10       -> Benign
    """

    def score_variant(
        self,
        gene: str,
        variant: str,
        evidence_codes: list[str] | None = None,
        evidence_bundle: EvidenceBundle | None = None,
    ) -> OncogenicityScoringResult:
        """Score a variant using explicit evidence codes or auto-detection.

        If evidence_codes is provided, uses those directly.
        If not, attempts auto-detection from the evidence bundle data.
        """
        auto_detected = False
        if evidence_codes:
            applied = self._apply_explicit_codes(evidence_codes)
        else:
            applied = self._auto_detect_codes(gene, variant, evidence_bundle)
            auto_detected = True

        total_points = sum(code.points for code in applied)
        classification = self._classify_by_points(total_points)
        confidence = self._assess_confidence(applied, auto_detected)

        return OncogenicityScoringResult(
            total_points=total_points,
            classification=classification,
            applied_codes=applied,
            auto_detected=auto_detected,
            confidence=confidence,
        )

    def explain_codes(self, codes: list[str]) -> str:
        """Return human-readable explanation of evidence codes."""
        lines = []
        for code in codes:
            info = ONCOGENICITY_EVIDENCE_CODES.get(code)
            if info:
                pts = int(info["points"])
                direction = "Oncogenic" if pts > 0 else "Benign"
                lines.append(
                    f"**{code}** ({info['strength']}, {pts:+d} points, "
                    f"{direction}): {info['description']}"
                )
            else:
                lines.append(f"**{code}**: Unknown evidence code")
        return "\n".join(lines)

    def _apply_explicit_codes(self, codes: list[str]) -> list[AppliedEvidenceCode]:
        """Apply explicitly provided evidence codes."""
        applied = []
        for code in codes:
            code = code.upper().strip()
            info = ONCOGENICITY_EVIDENCE_CODES.get(code)
            if info:
                applied.append(
                    AppliedEvidenceCode(
                        code=code,
                        points=int(info["points"]),
                        description=str(info["description"]),
                        evidence="Manually specified by user",
                    )
                )
            else:
                logger.warning("Unknown oncogenicity evidence code: %s", code)
        return applied

    def _auto_detect_codes(
        self,
        gene: str,
        variant: str,
        bundle: EvidenceBundle | None,
    ) -> list[AppliedEvidenceCode]:
        """Attempt to auto-detect applicable evidence codes from available data."""
        applied: list[AppliedEvidenceCode] = []
        gene_upper = gene.upper()

        # OVS1: Null variant in bona fide tumor suppressor gene
        if self._is_null_variant(variant) and gene_upper in KNOWN_TUMOR_SUPPRESSORS:
            applied.append(
                AppliedEvidenceCode(
                    code="OVS1",
                    points=8,
                    description="Null variant in bona fide tumor suppressor gene",
                    evidence=f"{variant} appears to be a null variant in {gene}, "
                    f"a known tumor suppressor",
                )
            )

        # OS3: Hotspot detection from CIViC/OncoKB
        if bundle and self._is_hotspot(gene, variant, bundle):
            applied.append(
                AppliedEvidenceCode(
                    code="OS3",
                    points=4,
                    description="Located in hotspot with sufficient statistical evidence",
                    evidence=f"Multiple evidence items support {gene} {variant} as a "
                    f"recurrent mutation hotspot",
                )
            )

        # OS1: Same amino acid change as established oncogenic variant
        if bundle and self._is_established_oncogenic(bundle):
            applied.append(
                AppliedEvidenceCode(
                    code="OS1",
                    points=4,
                    description="Same amino acid change as established oncogenic variant",
                    evidence="OncoKB or CIViC confirm this as an established oncogenic variant",
                )
            )

        # OS2: Functional studies from OncoKB
        if bundle and bundle.has_oncokb_data and bundle.oncokb_annotation:
            effect = bundle.oncokb_annotation.known_effect or ""
            if effect.lower() in (
                "gain-of-function",
                "loss-of-function",
                "switch-of-function",
            ):
                applied.append(
                    AppliedEvidenceCode(
                        code="OS2",
                        points=4,
                        description="Well-established functional studies show oncogenic effect",
                        evidence=f"OncoKB mutation effect: {effect}",
                    )
                )

        # OM1: Known functional domain (oncogene or TSG)
        if (gene_upper in KNOWN_ONCOGENES or gene_upper in KNOWN_TUMOR_SUPPRESSORS) and not any(
            c.code in ("OVS1", "OS3") for c in applied
        ):
            applied.append(
                AppliedEvidenceCode(
                    code="OM1",
                    points=2,
                    description="Located in critical/functional domain without benign variation",
                    evidence=f"{gene} is a known {'oncogene' if gene_upper in KNOWN_ONCOGENES else 'tumor suppressor'}",
                )
            )

        # OM4/OP4: Population frequency from ClinVar
        if bundle and bundle.has_clinvar_data:
            for cv in bundle.clinvar_variants:
                sig = (cv.clinical_significance or "").lower()
                if "benign" not in sig and not any(c.code in ("OM4", "OP4") for c in applied):
                    applied.append(
                        AppliedEvidenceCode(
                            code="OM4",
                            points=2,
                            description="Absent/extremely rare in population databases",
                            evidence="Not classified as benign in ClinVar, suggesting rare in population",
                        )
                    )

        # SBVS1/SBS1: Check for high population frequency (benign)
        if bundle and bundle.has_clinvar_data:
            for cv in bundle.clinvar_variants:
                sig = (cv.clinical_significance or "").lower()
                if "benign" in sig and "likely" not in sig:
                    # Strong benign signal
                    applied.append(
                        AppliedEvidenceCode(
                            code="SBS1",
                            points=-4,
                            description="MAF >1% in any general continental population in gnomAD",
                            evidence=f"ClinVar classifies as Benign (variation ID {cv.variation_id})",
                        )
                    )
                    break

        # OP2: Somatic variant in gene with single cancer etiology
        if (
            gene_upper in KNOWN_ONCOGENES
            and not any(c.code == "OP2" for c in applied)
            and len(applied) < 4
        ):
            applied.append(
                AppliedEvidenceCode(
                    code="OP2",
                    points=1,
                    description="Somatic variant in gene with single cancer etiology",
                    evidence=f"{gene} is a known oncogene",
                )
            )

        return applied

    @staticmethod
    def _is_null_variant(variant: str) -> bool:
        """Check if variant appears to be a null/truncating variant."""
        variant_lower = variant.lower().strip()
        # Check for known null variant patterns
        for null_type in NULL_VARIANT_TYPES:
            if null_type in variant_lower:
                return True
        # Check for frameshift notation (e.g., "fs", "Ter", "*")
        if re.search(r"(fs|ter\d*|\*\d*$|\*$|x\d+)", variant_lower):
            return True
        # Check for splice notation
        return bool(re.search(r"(splice|c\.\d+[+-]\d+)", variant_lower))

    @staticmethod
    def _is_hotspot(gene: str, variant: str, bundle: EvidenceBundle) -> bool:
        """Check if variant is in a known hotspot based on evidence volume."""
        # A variant with many CIViC evidence items likely represents a hotspot
        if len(bundle.civic_evidence) >= 5:
            return True
        # OncoKB oncogenic + high evidence level
        if bundle.has_oncokb_data and bundle.oncokb_annotation:
            oncogenic = (bundle.oncokb_annotation.oncogenic or "").lower()
            if "oncogenic" in oncogenic:
                return True
        return False

    @staticmethod
    def _is_established_oncogenic(bundle: EvidenceBundle) -> bool:
        """Check if variant is an established oncogenic variant."""
        if bundle.has_oncokb_data and bundle.oncokb_annotation:
            oncogenic = (bundle.oncokb_annotation.oncogenic or "").lower()
            if oncogenic == "oncogenic":
                return True
        # Check CIViC for multiple accepted oncogenic evidence items
        oncogenic_count = sum(
            1
            for e in bundle.civic_evidence
            if (e.significance or "").lower() in ("oncogenic", "gain of function")
        )
        return oncogenic_count >= 3

    @staticmethod
    def _classify_by_points(points: int) -> str:
        """Determine classification based on total points."""
        if points >= 10:
            return "Oncogenic"
        if points >= 6:
            return "Likely Oncogenic"
        if points <= -10:
            return "Benign"
        if points <= -6:
            return "Likely Benign"
        return "VUS"

    @staticmethod
    def _assess_confidence(applied: list[AppliedEvidenceCode], auto_detected: bool) -> str:
        """Assess confidence of the scoring."""
        if not applied:
            return "LOW"
        if auto_detected:
            if len(applied) >= 3:
                return "MODERATE"
            return "LOW"
        # Explicit codes: higher confidence
        if len(applied) >= 3:
            return "HIGH"
        if len(applied) >= 2:
            return "MODERATE"
        return "LOW"
