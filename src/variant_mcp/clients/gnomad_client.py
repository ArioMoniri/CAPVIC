"""gnomAD GraphQL client for population allele frequencies."""

from __future__ import annotations

import logging
from typing import Any

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.constants import GNOMAD_API_URL, GNOMAD_RATE_LIMIT
from variant_mcp.models.evidence import GnomADFrequency

logger = logging.getLogger(__name__)

# GraphQL query for gnomAD variant lookup
_VARIANT_QUERY = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    rsids
    flags
    genome {
      ac
      an
      af
      ac_hom
      filters
      populations {
        id
        ac
        an
        af
      }
    }
    exome {
      ac
      an
      af
      ac_hom
      filters
      populations {
        id
        ac
        an
        af
      }
    }
  }
}
"""

_SEARCH_QUERY = """
query GnomadSearch($query: String!, $datasetId: DatasetId!) {
  searchResults(query: $query, dataset: $datasetId) {
    label
    url
  }
}
"""

# Map genome version to gnomAD dataset ID
_DATASET_MAP = {
    "GRCh38": "gnomad_r4",
    "GRCh37": "gnomad_r2_1",
}

# Population IDs we care about
_POPULATION_IDS = {"afr", "amr", "asj", "eas", "fin", "nfe", "sas"}


class GnomADClient(BaseClient):
    """Client for gnomAD GraphQL API for population allele frequencies."""

    def __init__(self) -> None:
        super().__init__(
            base_url=GNOMAD_API_URL,
            rate_limit=GNOMAD_RATE_LIMIT,
            headers={"Content-Type": "application/json"},
        )

    async def get_variant_frequency(
        self, variant_id: str, genome_version: str = "GRCh38"
    ) -> GnomADFrequency:
        """Get population allele frequencies for a variant.

        Args:
            variant_id: Format like "1-55505647-T-C" (chrom-pos-ref-alt) or HGVS.
            genome_version: "GRCh38" (gnomAD v4) or "GRCh37" (gnomAD v2).

        Returns:
            GnomADFrequency with population allele frequency data.

        Raises:
            ClientError: If the variant is not found or the API returns an error.
        """
        dataset_id = _DATASET_MAP.get(genome_version)
        if dataset_id is None:
            raise ClientError(
                f"Unsupported genome version: {genome_version}. Use 'GRCh38' or 'GRCh37'."
            )

        payload = {
            "query": _VARIANT_QUERY,
            "variables": {"variantId": variant_id, "datasetId": dataset_id},
        }

        data = await self.post_json("", json_body=payload)
        errors = data.get("errors")
        if errors:
            msg = errors[0].get("message", str(errors))
            raise ClientError(f"gnomAD GraphQL error: {msg}")

        variant_data = (data.get("data") or {}).get("variant")
        if variant_data is None:
            raise ClientError(f"Variant {variant_id} not found in gnomAD ({genome_version}).")

        return self._parse_variant(variant_data, genome_version)

    async def search_by_gene_variant(self, gene: str, variant: str) -> GnomADFrequency | None:
        """Search by gene + protein change via the gnomAD search endpoint.

        Args:
            gene: Gene symbol, e.g. "BRAF".
            variant: Protein change, e.g. "V600E".

        Returns:
            GnomADFrequency if found, None otherwise.
        """
        query_str = f"{gene} {variant}"
        dataset_id = _DATASET_MAP["GRCh38"]

        payload = {
            "query": _SEARCH_QUERY,
            "variables": {"query": query_str, "datasetId": dataset_id},
        }

        data = await self.post_json("", json_body=payload)
        errors = data.get("errors")
        if errors:
            logger.warning("gnomAD search error: %s", errors[0].get("message"))
            return None

        results = (data.get("data") or {}).get("searchResults", [])
        if not results:
            return None

        # Extract variant ID from the first matching URL
        # gnomAD search returns URLs like "/variant/1-55505647-T-C?dataset=gnomad_r4"
        for result in results:
            url = result.get("url", "")
            if "/variant/" in url:
                vid = url.split("/variant/")[1].split("?")[0]
                try:
                    return await self.get_variant_frequency(vid)
                except ClientError:
                    continue

        return None

    @staticmethod
    def _parse_variant(data: dict[str, Any], genome_version: str) -> GnomADFrequency:
        """Parse gnomAD GraphQL variant response into a GnomADFrequency model."""
        # Prefer genome data, fall back to exome
        genome = data.get("genome") or {}
        exome = data.get("exome") or {}

        # Use genome AF if available, otherwise exome
        primary = genome if genome.get("af") is not None else exome

        af = primary.get("af")
        ac = primary.get("ac")
        an = primary.get("an")
        ac_hom = primary.get("ac_hom")

        # Combine population frequencies from genome and exome
        pop_freqs: dict[str, float] = {}
        for source in (genome, exome):
            for pop in source.get("populations", []):
                pop_id = pop.get("id", "")
                if (
                    pop_id in _POPULATION_IDS
                    and pop.get("af") is not None
                    and pop_id not in pop_freqs
                ):
                    pop_freqs[pop_id] = pop.get("af")

        # Collect filters
        filters = []
        for source in (genome, exome):
            src_filters = source.get("filters", [])
            if src_filters:
                filters.extend(src_filters)
        filter_status = ",".join(filters) if filters else "PASS"

        # Determine source label
        source = "gnomAD v4" if genome_version == "GRCh38" else "gnomAD v2"

        rsids = data.get("rsids", [])
        rsid = rsids[0] if rsids else None

        flags = data.get("flags", [])
        if flags:
            filter_status = f"{filter_status}; flags: {','.join(flags)}"

        return GnomADFrequency(
            variant_id=data.get("variant_id"),
            rsid=rsid,
            allele_frequency=af,
            allele_count=ac,
            allele_number=an,
            homozygote_count=ac_hom,
            population_frequencies=pop_freqs,
            genome_version=genome_version,
            filter_status=filter_status,
            source=source,
        )
