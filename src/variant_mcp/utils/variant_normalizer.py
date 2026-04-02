"""HGVS notation parser and variant normalizer — pure Python, no external dependencies."""

from __future__ import annotations

import re

from variant_mcp.models.evidence import VariantNotation


class VariantNormalizer:
    """Parse and normalize variant notation formats."""

    # Amino acid 1-letter to 3-letter mapping
    AA_MAP: dict[str, str] = {
        "A": "Ala",
        "R": "Arg",
        "N": "Asn",
        "D": "Asp",
        "C": "Cys",
        "E": "Glu",
        "Q": "Gln",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "L": "Leu",
        "K": "Lys",
        "M": "Met",
        "F": "Phe",
        "P": "Pro",
        "S": "Ser",
        "T": "Thr",
        "W": "Trp",
        "Y": "Tyr",
        "V": "Val",
        "*": "Ter",
        "X": "Ter",
    }

    # Reverse mapping: 3-letter to 1-letter
    AA_REVERSE: dict[str, str] = {v: k for k, v in AA_MAP.items() if k != "X"}

    # Pattern for 1-letter protein notation: V600E, p.V600E, K550fs, R175*
    _RE_PROTEIN_1LETTER = re.compile(
        r"^(?:p\.)?([A-Z*])(\d+)([A-Z*](?:fs\*?\d*)?|fs\*?\d*|del|dup|ins[A-Z]+)?$"
    )

    # Pattern for 3-letter protein notation: p.Val600Glu
    _RE_PROTEIN_3LETTER = re.compile(
        r"^(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}(?:fs[A-Z][a-z]{2}\*?\d*)?|del|dup)?$"
    )

    # Pattern for cDNA notation: c.1799T>A
    _RE_CDNA = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")

    # Pattern for cDNA indels: c.1234delA, c.1234_1235insAT, c.1234dupA
    _RE_CDNA_INDEL = re.compile(r"^c\.(\d+)(?:_(\d+))?(del[ACGT]*|dup[ACGT]*|ins[ACGT]+)$")

    # Pattern for splice site: c.1234+1G>A, c.1234-2A>C
    _RE_SPLICE = re.compile(r"^c\.(\d+)([+-]\d+)([ACGT])>([ACGT])$")

    def normalize(self, variant_input: str) -> VariantNotation:
        """Parse any variant notation format and return a structured VariantNotation."""
        cleaned = variant_input.strip()

        # Try 1-letter protein notation
        m = self._RE_PROTEIN_1LETTER.match(cleaned)
        if m:
            ref_aa, pos_str, alt = m.group(1), m.group(2), m.group(3)
            position = int(pos_str)
            alt_aa = alt if alt and len(alt) == 1 else alt
            variant_type = self._classify_protein_change(ref_aa, alt_aa)
            protein_1 = f"{ref_aa}{pos_str}{alt_aa}" if alt_aa else f"{ref_aa}{pos_str}"
            protein_3 = self._to_3letter_str(ref_aa, pos_str, alt_aa)
            return VariantNotation(
                original=variant_input,
                protein_1letter=protein_1,
                protein_3letter=f"p.{protein_3}",
                variant_type=variant_type,
                position=position,
                ref_aa=ref_aa,
                alt_aa=alt_aa if alt_aa and len(alt_aa) == 1 else None,
            )

        # Try 3-letter protein notation
        m = self._RE_PROTEIN_3LETTER.match(cleaned)
        if m:
            ref_3, pos_str, alt_3 = m.group(1), m.group(2), m.group(3)
            position = int(pos_str)
            ref_aa = self.AA_REVERSE.get(ref_3, "?")
            alt_aa = self.AA_REVERSE.get(alt_3, None) if alt_3 else None
            protein_1 = f"{ref_aa}{pos_str}{alt_aa}" if alt_aa else f"{ref_aa}{pos_str}"
            protein_3 = f"{ref_3}{pos_str}{alt_3}" if alt_3 else f"{ref_3}{pos_str}"
            variant_type = self._classify_protein_change(ref_aa, alt_aa)
            return VariantNotation(
                original=variant_input,
                protein_1letter=protein_1,
                protein_3letter=f"p.{protein_3}",
                variant_type=variant_type,
                position=position,
                ref_aa=ref_aa,
                alt_aa=alt_aa,
            )

        # Try cDNA substitution
        m = self._RE_CDNA.match(cleaned)
        if m:
            pos_str = m.group(1)
            return VariantNotation(
                original=variant_input,
                cdna=cleaned,
                variant_type="substitution",
                position=int(pos_str),
            )

        # Try cDNA indel
        m = self._RE_CDNA_INDEL.match(cleaned)
        if m:
            pos_str = m.group(1)
            change = m.group(3)
            if change.startswith("del"):
                vtype = "deletion"
            elif change.startswith("dup"):
                vtype = "duplication"
            else:
                vtype = "insertion"
            return VariantNotation(
                original=variant_input,
                cdna=cleaned,
                variant_type=vtype,
                position=int(pos_str),
            )

        # Try splice site
        m = self._RE_SPLICE.match(cleaned)
        if m:
            pos_str = m.group(1)
            return VariantNotation(
                original=variant_input,
                cdna=cleaned,
                variant_type="splice",
                position=int(pos_str),
            )

        # Fallback: return what we can
        fallback_pos = self.extract_position(cleaned)
        return VariantNotation(
            original=variant_input,
            position=fallback_pos,
        )

    def to_protein_1letter(self, variant: str) -> str:
        """Convert to 1-letter protein change: V600E"""
        parsed = self.normalize(variant)
        if parsed.protein_1letter:
            return parsed.protein_1letter
        return variant

    def to_protein_3letter(self, variant: str) -> str:
        """Convert to 3-letter HGVS protein: p.Val600Glu"""
        parsed = self.normalize(variant)
        if parsed.protein_3letter:
            return parsed.protein_3letter
        return variant

    def to_short(self, variant: str) -> str:
        """Convert to the shortest standard form: V600E"""
        return self.to_protein_1letter(variant)

    def detect_variant_type(self, variant: str) -> str:
        """Detect type: missense, nonsense, frameshift, splice, indel, etc."""
        parsed = self.normalize(variant)
        return parsed.variant_type or "unknown"

    def extract_position(self, variant: str) -> int | None:
        """Extract amino acid position number from the variant notation."""
        m = re.search(r"(\d+)", variant)
        if m:
            return int(m.group(1))
        return None

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _classify_protein_change(self, ref_aa: str | None, alt_aa: str | None) -> str:
        """Classify a protein-level change based on reference and alternate amino acids."""
        if alt_aa is None:
            return "unknown"
        if isinstance(alt_aa, str) and ("fs" in alt_aa.lower()):
            return "frameshift"
        if alt_aa in ("del",):
            return "deletion"
        if alt_aa in ("dup",):
            return "duplication"
        if isinstance(alt_aa, str) and alt_aa.startswith("ins"):
            return "insertion"
        if alt_aa in ("*", "X") or alt_aa == "Ter":
            return "nonsense"
        if ref_aa == alt_aa:
            return "synonymous"
        if ref_aa and alt_aa and len(alt_aa) == 1:
            return "missense"
        return "unknown"

    def _to_3letter_str(self, ref_1: str, pos: str, alt_1: str | None) -> str:
        """Convert 1-letter ref/alt to 3-letter string."""
        ref_3 = self.AA_MAP.get(ref_1, ref_1)
        if alt_1 is None:
            return f"{ref_3}{pos}"
        if alt_1 in ("del", "dup") or alt_1.startswith("ins") or "fs" in alt_1:
            return f"{ref_3}{pos}{alt_1}"
        alt_3 = self.AA_MAP.get(alt_1, alt_1)
        return f"{ref_3}{pos}{alt_3}"
