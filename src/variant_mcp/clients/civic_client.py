"""CIViC V2 GraphQL API client."""

from __future__ import annotations

import logging
from typing import Any

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.constants import CIVIC_BASE_URL, CIVIC_GRAPHQL_URL, CIVIC_RATE_LIMIT
from variant_mcp.models.evidence import CIViCAssertion, CIViCEvidenceItem, SourceInfo
from variant_mcp.queries.civic_graphql import CIViCQueries

logger = logging.getLogger(__name__)


class CIViCClient(BaseClient):
    """Client for the CIViC V2 GraphQL API."""

    def __init__(self) -> None:
        super().__init__(
            base_url=CIVIC_GRAPHQL_URL,
            rate_limit=CIVIC_RATE_LIMIT,
            headers={"Content-Type": "application/json"},
        )

    async def _graphql(self, query: str, variables: dict[str, Any] | None = None) -> dict:
        """Execute a GraphQL query."""
        body: dict[str, Any] = {"query": query}
        if variables:
            body["variables"] = {k: v for k, v in variables.items() if v is not None}
        result = await self.post_json("", json_body=body)
        if "errors" in result:
            raise ClientError(f"CIViC GraphQL error: {result['errors']}")
        return result.get("data", {})

    async def search_evidence(
        self,
        *,
        disease: str | None = None,
        gene: str | None = None,
        variant: str | None = None,
        therapy: str | None = None,
        evidence_type: str | None = None,
        significance: str | None = None,
        first: int = 25,
        after: str | None = None,
    ) -> dict[str, Any]:
        """Search CIViC evidence items with filters."""
        variables = {
            "diseaseName": disease,
            "geneName": gene,
            "variantName": variant,
            "therapyName": therapy,
            "evidenceType": evidence_type,
            "significance": significance,
            "first": first,
            "after": after,
        }
        data = await self._graphql(CIViCQueries.SEARCH_EVIDENCE, variables)
        return data.get("evidenceItems", {})

    async def search_evidence_parsed(
        self,
        *,
        disease: str | None = None,
        gene: str | None = None,
        variant: str | None = None,
        therapy: str | None = None,
        evidence_type: str | None = None,
        first: int = 25,
    ) -> list[CIViCEvidenceItem]:
        """Search and return parsed CIViC evidence items."""
        raw = await self.search_evidence(
            disease=disease, gene=gene, variant=variant,
            therapy=therapy, evidence_type=evidence_type, first=first,
        )
        items = []
        for node in raw.get("nodes", []):
            items.append(self._parse_evidence_node(node))
        return items

    async def get_gene(self, name: str) -> dict[str, Any]:
        """Get gene details and associated variants."""
        data = await self._graphql(CIViCQueries.GET_GENE, {"name": name})
        genes = data.get("genes", {}).get("nodes", [])
        return genes[0] if genes else {}

    async def get_variant(self, variant_id: int) -> dict[str, Any]:
        """Get variant details by CIViC ID."""
        data = await self._graphql(CIViCQueries.GET_VARIANT, {"id": variant_id})
        return data.get("variant", {})

    async def get_evidence_item(self, evidence_id: int) -> CIViCEvidenceItem | None:
        """Get a single evidence item by ID."""
        data = await self._graphql(CIViCQueries.GET_EVIDENCE_ITEM, {"id": evidence_id})
        node = data.get("evidenceItem")
        if node:
            return self._parse_evidence_node(node)
        return None

    async def search_assertions(
        self,
        *,
        disease: str | None = None,
        gene: str | None = None,
        therapy: str | None = None,
        significance: str | None = None,
        first: int = 25,
        after: str | None = None,
    ) -> list[CIViCAssertion]:
        """Search CIViC curated assertions."""
        variables = {
            "diseaseName": disease,
            "geneName": gene,
            "therapyName": therapy,
            "significance": significance,
            "first": first,
            "after": after,
        }
        data = await self._graphql(CIViCQueries.SEARCH_ASSERTIONS, variables)
        assertions_data = data.get("assertions", {})
        results = []
        for node in assertions_data.get("nodes", []):
            results.append(self._parse_assertion_node(node))
        return results

    async def search_typeahead(
        self, entity_type: str, query_term: str
    ) -> list[dict[str, Any]]:
        """Autocomplete search for genes, diseases, or therapies."""
        query_map = {
            "gene": CIViCQueries.TYPEAHEAD_GENES,
            "disease": CIViCQueries.TYPEAHEAD_DISEASES,
            "therapy": CIViCQueries.TYPEAHEAD_THERAPIES,
        }
        query = query_map.get(entity_type.lower())
        if not query:
            raise ClientError(
                f"Unknown entity type '{entity_type}'. Use 'gene', 'disease', or 'therapy'."
            )
        data = await self._graphql(query, {"queryTerm": query_term})
        key = f"{entity_type.lower()}Typeahead"
        return data.get(key, [])

    @staticmethod
    def _parse_evidence_node(node: dict) -> CIViCEvidenceItem:
        mp = node.get("molecularProfile", {}) or {}
        variants = mp.get("variants", []) or []
        disease_data = node.get("disease", {}) or {}
        therapies = node.get("therapies", []) or []
        phenotypes = node.get("phenotypes", []) or []
        source_data = node.get("source", {}) or {}

        source = SourceInfo(
            source_type=source_data.get("sourceType"),
            citation=source_data.get("citation"),
            source_url=source_data.get("sourceUrl"),
            pmid=source_data.get("citationId"),
        ) if source_data else None

        return CIViCEvidenceItem(
            id=node["id"],
            name=node.get("name"),
            status=node.get("status"),
            evidence_type=node.get("evidenceType"),
            evidence_level=node.get("evidenceLevel"),
            evidence_direction=node.get("evidenceDirection"),
            evidence_rating=node.get("evidenceRating"),
            significance=node.get("significance"),
            description=node.get("description"),
            therapy_interaction_type=node.get("therapyInteractionType"),
            gene=variants[0].get("name", "").split()[0] if variants else None,
            variant=variants[0].get("name") if variants else None,
            molecular_profile=mp.get("name"),
            disease=disease_data.get("name"),
            disease_doid=disease_data.get("doid"),
            therapies=[t.get("name", "") for t in therapies],
            phenotypes=[p.get("name", "") for p in phenotypes],
            source=source,
        )

    @staticmethod
    def _parse_assertion_node(node: dict) -> CIViCAssertion:
        mp = node.get("molecularProfile", {}) or {}
        variants = mp.get("variants", []) or []
        disease_data = node.get("disease", {}) or {}
        therapies = node.get("therapies", []) or []
        nccn = node.get("nccnGuideline", {}) or {}
        acmg = node.get("acmgCodes", []) or []
        evidence_items = node.get("evidenceItems", []) or []

        return CIViCAssertion(
            id=node["id"],
            name=node.get("name"),
            assertion_type=node.get("assertionType"),
            assertion_direction=node.get("assertionDirection"),
            significance=node.get("significance"),
            summary=node.get("summary"),
            description=node.get("description"),
            amp_level=node.get("ampLevel"),
            gene=variants[0].get("name", "").split()[0] if variants else None,
            variant=variants[0].get("name") if variants else None,
            molecular_profile=mp.get("name"),
            disease=disease_data.get("name"),
            therapies=[t.get("name", "") for t in therapies],
            nccn_guideline=nccn.get("name"),
            acmg_codes=[c.get("code", "") for c in acmg],
            evidence_count=len(evidence_items),
        )
