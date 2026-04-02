"""Summary table formatters."""

from __future__ import annotations

from variant_mcp.constants import (
    ACMG_CRITERIA,
    AMP_TIER_DEFINITIONS,
    CIVIC_EVIDENCE_LEVELS,
    DISCLAIMER,
    ONCOGENICITY_EVIDENCE_CODES,
    ONCOKB_LEVELS,
    REFERENCE_PMIDS,
)


class TableFormatter:
    """Generates formatted reference tables for classification frameworks."""

    @staticmethod
    def format_frameworks_reference() -> str:
        """Comprehensive reference document for all classification frameworks."""
        lines: list[str] = []
        lines.append("# Classification Frameworks Reference\n")

        # AMP/ASCO/CAP
        lines.append("## 1. AMP/ASCO/CAP 4-Tier Somatic Classification\n")
        lines.append(f"*Li et al., J Mol Diagn 2017. PMID: {REFERENCE_PMIDS['AMP_ASCO_CAP']}*\n")
        lines.append("**Use for**: Somatic (cancer) variant clinical significance\n")
        lines.append("| Tier | Name | Evidence Level | Description |")
        lines.append("|------|------|----------------|-------------|")
        for tier, info in AMP_TIER_DEFINITIONS.items():
            name = info["name"]
            if "level_A" in info:
                lines.append(f"| {tier} | {name} | Level A | {info['level_A']} |")
                lines.append(f"| {tier} | {name} | Level B | {info['level_B']} |")
            elif "level_C" in info:
                lines.append(f"| {tier} | {name} | Level C | {info['level_C']} |")
                lines.append(f"| {tier} | {name} | Level D | {info['level_D']} |")
            else:
                lines.append(f"| {tier} | {name} | — | {info.get('description', '')} |")
        lines.append("")

        # ClinGen/CGC/VICC Oncogenicity SOP
        lines.append("## 2. ClinGen/CGC/VICC Oncogenicity SOP\n")
        lines.append(
            f"*Horak et al., Genet Med 2022. PMID: {REFERENCE_PMIDS['ONCOGENICITY_SOP']}*\n"
        )
        lines.append("**Use for**: Somatic variant oncogenicity assessment (point-based)\n")
        lines.append("### Evidence Codes\n")
        lines.append("| Code | Strength | Points | Description |")
        lines.append("|------|----------|--------|-------------|")
        for code, onco_info in ONCOGENICITY_EVIDENCE_CODES.items():
            points = int(onco_info["points"])
            lines.append(
                f"| {code} | {onco_info['strength']} | {points:+d} | {onco_info['description']} |"
            )
        lines.append("\n### Classification Thresholds\n")
        lines.append("| Classification | Points |")
        lines.append("|---------------|--------|")
        lines.append("| Oncogenic | >= 10 |")
        lines.append("| Likely Oncogenic | 6-9 |")
        lines.append("| VUS | -5 to 5 |")
        lines.append("| Likely Benign | -6 to -9 |")
        lines.append("| Benign | <= -10 |")
        lines.append("")

        # ACMG/AMP Germline
        lines.append("## 3. ACMG/AMP 5-Tier Germline Classification\n")
        lines.append(f"*Richards et al., Genet Med 2015. PMID: {REFERENCE_PMIDS['ACMG_AMP']}*\n")
        lines.append("**Use for**: Germline variant pathogenicity\n")
        lines.append(
            "**Classifications**: Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign\n"
        )
        lines.append("### Key Criteria\n")
        lines.append("| Code | Description |")
        lines.append("|------|-------------|")
        for code, desc in ACMG_CRITERIA.items():
            lines.append(f"| {code} | {desc} |")
        lines.append("")

        # When to use which framework
        lines.append("## When to Use Which Framework\n")
        lines.append("| Context | Framework | Notes |")
        lines.append("|---------|-----------|-------|")
        lines.append(
            "| Somatic variant, clinical actionability | AMP/ASCO/CAP | "
            "Focus on therapeutic, diagnostic, prognostic significance |"
        )
        lines.append(
            "| Somatic variant, oncogenicity question | ClinGen/CGC/VICC SOP | "
            "Point-based scoring for oncogenic vs benign |"
        )
        lines.append(
            "| Germline variant, disease causation | ACMG/AMP | "
            "For inherited variants, requires expert review |"
        )
        lines.append(
            "| Both somatic frameworks | AMP + Oncogenicity SOP | "
            "These are complementary — AMP focuses on clinical action, SOP on biology |"
        )

        # Evidence sources
        lines.append("\n## Data Sources\n")
        lines.append("| Source | Data Type | Access |")
        lines.append("|--------|-----------|--------|")
        lines.append("| CIViC | Expert-curated clinical evidence | Free, no auth |")
        lines.append("| ClinVar | Aggregate pathogenicity from labs | Free, no auth |")
        lines.append("| OncoKB | FDA-recognized oncogenicity | Free academic token |")
        lines.append("| VICC MetaKB | Harmonized from 6 KBs | Free, no auth |")

        lines.append(f"\n---\n\n{DISCLAIMER}")
        return "\n".join(lines)

    @staticmethod
    def format_evidence_levels_table() -> str:
        """Format CIViC evidence levels reference table."""
        lines = ["## CIViC Evidence Levels\n"]
        lines.append("| Level | Description |")
        lines.append("|-------|-------------|")
        for level, desc in CIVIC_EVIDENCE_LEVELS.items():
            lines.append(f"| {level} | {desc} |")
        return "\n".join(lines)

    @staticmethod
    def format_oncokb_levels_table() -> str:
        """Format OncoKB therapeutic levels reference table."""
        lines = ["## OncoKB Therapeutic Levels of Evidence\n"]
        lines.append("| Level | Description |")
        lines.append("|-------|-------------|")
        for level, desc in ONCOKB_LEVELS.items():
            lines.append(f"| {level} | {desc} |")
        return "\n".join(lines)
