"""VICC MetaKB search client."""

from __future__ import annotations

import logging
from typing import Any

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.constants import METAKB_RATE_LIMIT, METAKB_SEARCH_URL
from variant_mcp.models.evidence import MetaKBInterpretation

logger = logging.getLogger(__name__)

METAKB_FALLBACK_MSG = (
    "VICC MetaKB REST API was not reachable. You can manually cross-reference "
    "results at https://search.cancervariants.org. Results from CIViC and ClinVar "
    "are still available and cover most MetaKB member knowledgebases."
)


class MetaKBClient(BaseClient):
    """Client for the VICC MetaKB search API.

    Falls back gracefully if the REST API endpoint is not accessible.
    """

    def __init__(self) -> None:
        super().__init__(
            base_url=METAKB_SEARCH_URL,
            rate_limit=METAKB_RATE_LIMIT,
            headers={"Accept": "application/json"},
        )

    async def search(
        self,
        *,
        gene: str | None = None,
        variant: str | None = None,
        disease: str | None = None,
    ) -> list[MetaKBInterpretation]:
        """Search MetaKB for harmonized interpretations.

        Returns parsed interpretations or empty list with logged warning on failure.
        """
        query_parts = []
        if gene:
            query_parts.append(gene)
        if variant:
            query_parts.append(variant)
        if disease:
            query_parts.append(disease)

        if not query_parts:
            return []

        try:
            data = await self._search_api(query_parts)
            return self._parse_results(data)
        except (ClientError, Exception) as e:
            logger.warning("MetaKB search failed (falling back): %s", e)
            return []

    async def _search_api(self, query_parts: list[str]) -> dict[str, Any]:
        """Attempt MetaKB API search."""
        # The MetaKB API endpoint structure — try the /search endpoint
        search_term = " ".join(query_parts)
        try:
            data = await self.get_json("api/v2/search", params={"q": search_term})
            return data  # type: ignore[no-any-return]
        except ClientError:
            # Try alternative endpoint format
            try:
                data = await self.get_json("api/search", params={"term": search_term})
                return data  # type: ignore[no-any-return]
            except ClientError as exc:
                raise ClientError(METAKB_FALLBACK_MSG) from exc

    @staticmethod
    def _parse_results(data: dict[str, Any]) -> list[MetaKBInterpretation]:
        """Parse MetaKB search results into MetaKBInterpretation list."""
        interpretations = []

        # Handle various response formats from MetaKB
        results = data.get("results", data.get("matches", data.get("data", [])))
        if not isinstance(results, list):
            results = []

        for result in results:
            source = result.get("source", result.get("src", ""))
            disease = result.get("disease", result.get("condition", ""))
            if isinstance(disease, dict):
                disease = disease.get("name", str(disease))

            drugs = result.get("drugs", result.get("therapies", []))
            if isinstance(drugs, list):
                drugs = [d.get("name", str(d)) if isinstance(d, dict) else str(d) for d in drugs]
            else:
                drugs = []

            interpretations.append(
                MetaKBInterpretation(
                    source=source,
                    disease=disease,
                    drugs=drugs,
                    evidence_level=result.get("evidence_level", result.get("tier")),
                    clinical_significance=result.get(
                        "clinical_significance",
                        result.get("significance"),
                    ),
                    description=result.get("description", result.get("summary")),
                    url=result.get("url", result.get("link")),
                )
            )
        return interpretations
