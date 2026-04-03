"""LitVar2 API client (NCBI).

Provides variant-centric literature mining — maps variants to publications,
co-mentioned diseases, clinical significance, and HGVS nomenclatures.

References:
  - LitVar2: Allot et al. 2023, Nucleic Acids Res. PMID: 36350631
  - API: https://www.ncbi.nlm.nih.gov/research/litvar2/
"""

from __future__ import annotations

import logging
from typing import Any

from variant_mcp.clients.base_client import BaseClient
from variant_mcp.constants import LITVAR2_API_URL, LITVAR_RATE_LIMIT, LITVAR_V1_API_URL

logger = logging.getLogger(__name__)


class LitVarClient(BaseClient):
    """Client for the NCBI LitVar2 API."""

    def __init__(self) -> None:
        super().__init__(
            base_url=LITVAR2_API_URL,
            rate_limit=LITVAR_RATE_LIMIT,
            headers={"Accept": "application/json"},
        )

    async def autocomplete(self, query: str) -> list[dict[str, Any]]:
        """Search LitVar2 for variants matching a gene/variant query.

        Args:
            query: Free-text query (e.g., "BRAF V600E", "TP53 R175H", "rs113488022").

        Returns:
            List of matching variant records with rsid, gene, pmids_count, etc.
        """
        try:
            result = await self.get_json("variant/autocomplete/", params={"query": query})
            return result if isinstance(result, list) else []
        except Exception as e:
            logger.warning("LitVar2 autocomplete failed for %r: %s", query, e)
            return []

    async def get_variant(self, litvar_id: str) -> dict[str, Any] | None:
        """Get full variant details from LitVar2.

        Args:
            litvar_id: LitVar ID (e.g., "litvar@rs113488022##").

        Returns:
            Variant detail dict or None if not found.
        """
        try:
            # URL-encode '#' chars — they are URL fragment separators and httpx strips them
            # No trailing slash — LitVar2 returns 404 with trailing slash on this endpoint
            safe_id = litvar_id.replace("#", "%23")
            result = await self.get_json(f"variant/get/{safe_id}")
            return result if isinstance(result, dict) else None
        except Exception as e:
            logger.warning("LitVar2 variant get failed for %s: %s", litvar_id, e)
            return None

    async def get_variant_entity(self, rsid: str) -> dict[str, Any] | None:
        """Get enriched variant entity from LitVar v1 API.

        The v1 entity endpoint returns richer data: disease co-mentions with
        counts, year distribution, all HGVS forms with frequencies, and gene
        details. Falls back gracefully if unavailable.

        Args:
            rsid: dbSNP rsID (e.g., "rs113488022"). Do NOT include "rs" prefix
                  if already in the ID.

        Returns:
            Enriched variant entity dict or None.
        """
        clean_rsid = rsid.lstrip("rs") if rsid.startswith("rs") else rsid
        url = f"{LITVAR_V1_API_URL}/entity/litvar/rs{clean_rsid}%23%23"
        try:
            client = await self._get_client()
            await self._rate_limiter.acquire()
            response = await client.get(url)
            response.raise_for_status()
            data = response.json()
            return data if isinstance(data, dict) else None
        except Exception as e:
            logger.warning("LitVar v1 entity failed for rs%s: %s", clean_rsid, e)
            return None

    async def search_variant_literature(self, gene: str, variant: str) -> dict[str, Any]:
        """High-level search: find variant in LitVar2, then enrich with v1 entity.

        Returns a combined result with publication count, diseases, clinical
        significance, and HGVS forms.
        """
        # Step 1: Autocomplete search
        matches = await self.autocomplete(f"{gene} {variant}")
        if not matches:
            return {
                "found": False,
                "gene": gene,
                "variant": variant,
                "pmids_count": 0,
            }

        best = matches[0]
        result: dict[str, Any] = {
            "found": True,
            "gene": gene,
            "variant": variant,
            "litvar_id": best.get("_id"),
            "rsid": best.get("rsid"),
            "name": best.get("name"),
            "hgvs": best.get("hgvs"),
            "pmids_count": best.get("pmids_count", 0),
            "clinical_significance": best.get("data_clinical_significance", []),
            "genes": best.get("gene", []),
        }

        # Step 2: Get LitVar2 full detail
        litvar_id = best.get("_id")
        if litvar_id:
            detail = await self.get_variant(litvar_id)
            if detail:
                result["clingen_ids"] = detail.get("clingen_ids", [])
                result["chromosome_position"] = detail.get("data_chromosome_base_position", [])
                result["snp_class"] = detail.get("data_snp_class", [])

        # Step 3: Enrich with v1 entity (disease co-mentions, HGVS forms)
        rsid = best.get("rsid")
        if rsid:
            entity = await self.get_variant_entity(rsid)
            if entity:
                result["diseases"] = entity.get("diseases", {})
                result["years"] = entity.get("years", [])
                result["all_hgvs"] = [
                    {"notation": h[0], "count": h[1]}
                    for h in entity.get("all_hgvs", [])
                    if isinstance(h, list) and len(h) == 2
                ]
                result["first_published_year"] = entity.get("first_published_year")
                # Override pmids_count with v1 count if higher (v1 may have more)
                v1_count = entity.get("pmids_count", 0)
                if v1_count > result.get("pmids_count", 0):
                    result["pmids_count"] = v1_count

        return result
