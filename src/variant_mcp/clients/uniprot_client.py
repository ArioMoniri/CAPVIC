"""UniProt REST API client for protein domain and feature information."""

from __future__ import annotations

import logging
import sys
from typing import Any

from variant_mcp.clients.base_client import BaseClient
from variant_mcp.constants import UNIPROT_API_URL, UNIPROT_RATE_LIMIT
from variant_mcp.models.evidence import (
    DomainCheckResult,
    ProteinDomain,
    ProteinFeatures,
)
from variant_mcp.utils.variant_normalizer import VariantNormalizer

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stderr))


class UniProtClient(BaseClient):
    """Client for UniProt REST API for protein domain information."""

    def __init__(self) -> None:
        super().__init__(
            base_url=UNIPROT_API_URL,
            rate_limit=UNIPROT_RATE_LIMIT,
            headers={"Accept": "application/json"},
        )
        self._normalizer = VariantNormalizer()

    async def get_protein_features(self, gene: str) -> ProteinFeatures:
        """Get protein features (domains, active sites, etc.) for a gene.

        Only retrieves human (organism_id=9606) proteins.
        """
        params = {
            "query": f"gene_exact:{gene} AND organism_id:9606",
            "fields": (
                "ft_domain,ft_region,ft_act_site,ft_binding,ft_site,"
                "ft_mutagen,length,accession,gene_names"
            ),
            "format": "json",
            "size": "1",
        }
        data = await self.get_json("uniprotkb/search", params=params)
        results: list[dict[str, Any]] = data.get("results", [])

        if not results:
            return ProteinFeatures(gene=gene)

        entry = results[0]
        uniprot_id: str | None = entry.get("primaryAccession")
        length: int | None = entry.get("sequence", {}).get("length")
        domains = self._extract_domains(entry)

        return ProteinFeatures(
            gene=gene,
            uniprot_id=uniprot_id,
            protein_length=length,
            domains=domains,
        )

    async def get_domain_at_position(self, gene: str, position: int) -> list[ProteinDomain]:
        """Check which functional domains overlap a given amino acid position."""
        features = await self.get_protein_features(gene)
        return [d for d in features.domains if d.start_pos <= position <= d.end_pos]

    async def check_variant_in_domain(self, gene: str, variant: str) -> DomainCheckResult:
        """Given gene + variant, check if the variant position falls in a known functional domain.

        This is the key method for OM1 evidence code in oncogenicity SOP.
        """
        parsed = self._normalizer.normalize(variant)
        position = parsed.position

        if position is None:
            return DomainCheckResult(gene=gene, variant=variant)

        overlapping = await self.get_domain_at_position(gene, position)

        return DomainCheckResult(
            gene=gene,
            variant=variant,
            position=position,
            in_domain=bool(overlapping),
            domains=overlapping,
            evidence_for_om1=bool(overlapping),
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_domains(entry: dict[str, Any]) -> list[ProteinDomain]:
        """Extract domain/feature annotations from a UniProt JSON entry."""
        domains: list[ProteinDomain] = []
        features: list[dict[str, Any]] = entry.get("features", [])

        for feat in features:
            feat_type = feat.get("type", "")
            location = feat.get("location", {})
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")

            if start is None or end is None:
                continue

            description = feat.get("description", "")
            source = None
            evidences = feat.get("evidences", [])
            for ev in evidences:
                src = ev.get("source", {}).get("name")
                if src:
                    source = src
                    break

            domains.append(
                ProteinDomain(
                    name=description or feat_type,
                    domain_type=feat_type,
                    start_pos=int(start),
                    end_pos=int(end),
                    source=source,
                    description=description or None,
                )
            )

        return domains
