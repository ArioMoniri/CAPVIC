"""Cancer Hotspots API client (cancerhotspots.org).

Provides recurrent somatic mutation hotspot data from large-scale cancer
genomics studies (Chang et al. 2016, 2018).
"""

from __future__ import annotations

import logging
from typing import Any

from variant_mcp.clients.base_client import BaseClient
from variant_mcp.constants import CANCER_HOTSPOTS_API_URL, CANCER_HOTSPOTS_RATE_LIMIT

logger = logging.getLogger(__name__)


class CancerHotspotsClient(BaseClient):
    """Client for the Cancer Hotspots API (cancerhotspots.org)."""

    def __init__(self) -> None:
        super().__init__(
            base_url=CANCER_HOTSPOTS_API_URL,
            rate_limit=CANCER_HOTSPOTS_RATE_LIMIT,
            headers={"Accept": "application/json"},
        )

    async def get_hotspots_by_gene(self, gene: str) -> list[dict[str, Any]]:
        """Get all hotspot residues for a gene.

        The Cancer Hotspots API requires POST with a JSON array body.
        Returns list of hotspot residues with sample counts per cancer type.
        """
        try:
            response = await self.post("hotspots/single/byGene", json_body=[gene])  # type: ignore[arg-type]
            result = response.json()
            return result if isinstance(result, list) else []
        except Exception as e:
            logger.warning("Cancer Hotspots gene lookup failed for %s: %s", gene, e)
            return []

    async def get_hotspot_at_residue(self, gene: str, residue: int) -> list[dict[str, Any]]:
        """Check if a specific residue is a known hotspot."""
        hotspots = await self.get_hotspots_by_gene(gene)
        return [h for h in hotspots if h.get("residue") and str(residue) in str(h["residue"])]

    async def get_3d_hotspots_by_gene(self, gene: str) -> list[dict[str, Any]]:
        """Get 3D structure-based hotspot clusters for a gene.

        Returns clusters identified by 3D proximity in protein structure.
        """
        try:
            response = await self.post("hotspots/3d/byGene", json_body=[gene])  # type: ignore[arg-type]
            result = response.json()
            return result if isinstance(result, list) else []
        except Exception as e:
            logger.warning("Cancer Hotspots 3D lookup failed for %s: %s", gene, e)
            return []

    def parse_hotspot_residue(self, hotspot: dict) -> dict[str, Any]:
        """Parse a hotspot record into a standardized format."""
        residue = hotspot.get("residue", "")
        variant_counts = hotspot.get("variantAminoAcid", {}) or {}
        total_count = sum(int(v) for v in variant_counts.values() if str(v).isdigit())
        cancer_types = hotspot.get("tumorTypeComposition", {}) or {}
        # Parse q-value — API returns string or None
        raw_q = hotspot.get("qValue")
        q_value = float(raw_q) if raw_q is not None else None
        raw_qct = hotspot.get("qValueCancerType")
        q_value_ct = float(raw_qct) if raw_qct is not None else None

        return {
            "gene": hotspot.get("hugoSymbol"),
            "residue": residue,
            "total_sample_count": total_count,
            "variant_amino_acids": variant_counts,
            "cancer_type_counts": cancer_types,
            "classification": hotspot.get("type", ""),
            "q_value": q_value,
            "q_value_cancertype": q_value_ct,
        }
