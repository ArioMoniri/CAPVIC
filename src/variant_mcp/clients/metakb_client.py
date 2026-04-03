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
        except (ClientError, KeyError, TypeError, ValueError) as e:
            logger.warning("MetaKB search failed (falling back): %s", e)
            return []

    async def _search_api(self, query_parts: list[str]) -> dict[str, Any]:
        """Attempt MetaKB API search.

        MetaKB v1 API endpoint: /api/v1/associations?q=<term>&size=<n>
        Swagger: https://search.cancervariants.org/api/v1/ui/
        """
        search_term = " ".join(query_parts)
        try:
            data = await self.get_json(
                "api/v1/associations", params={"q": search_term, "size": "25"}
            )
            return data  # type: ignore[no-any-return]
        except ClientError as exc:
            raise ClientError(METAKB_FALLBACK_MSG) from exc

    @staticmethod
    def _parse_results(data: dict[str, Any]) -> list[MetaKBInterpretation]:
        """Parse MetaKB v1 search results into MetaKBInterpretation list.

        MetaKB v1 /associations returns: {"hits": {"hits": [{"association": {...}}]}}
        Each association has: description, evidence_label, evidence_level,
        phenotypes, environmentalContexts, response_type, source_link, variant_name.
        """
        interpretations = []

        # MetaKB v1 response: {"hits": {"hits": [{"association": {...}}]}}
        hits_outer = data.get("hits", {})
        hits = hits_outer.get("hits", []) if isinstance(hits_outer, dict) else []

        for hit in hits:
            assoc = hit.get("association", hit) if isinstance(hit, dict) else {}

            # Extract disease from phenotypes list
            phenotypes = assoc.get("phenotypes", [])
            disease = ""
            if phenotypes and isinstance(phenotypes, list):
                disease = (
                    phenotypes[0].get("description", "") if isinstance(phenotypes[0], dict) else ""
                )
            if not disease:
                disease = assoc.get("disease_labels_truncated", "")

            # Extract drugs from environmentalContexts
            env_contexts = assoc.get("environmentalContexts", [])
            drugs = []
            if isinstance(env_contexts, list):
                drugs = [
                    ctx.get("description", str(ctx)) if isinstance(ctx, dict) else str(ctx)
                    for ctx in env_contexts
                    if ctx
                ]

            # Extract source from evidence list or source_link
            source = ""
            evidence = assoc.get("evidence", [])
            if evidence and isinstance(evidence, list) and isinstance(evidence[0], dict):
                ev_type = evidence[0].get("evidenceType", {})
                source = ev_type.get("sourceName", "") if isinstance(ev_type, dict) else ""

            # Map evidence_label to level (A-E style from CIViC)
            evidence_label = assoc.get("evidence_label", "")

            interpretations.append(
                MetaKBInterpretation(
                    source=source,
                    disease=disease,
                    drugs=drugs,
                    evidence_level=evidence_label or assoc.get("evidence_level"),
                    clinical_significance=assoc.get("response_type", ""),
                    description=assoc.get("description", ""),
                    url=assoc.get("source_link", ""),
                )
            )
        return interpretations
