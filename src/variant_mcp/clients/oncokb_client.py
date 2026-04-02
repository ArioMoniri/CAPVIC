"""OncoKB REST API client."""

from __future__ import annotations

import logging
import os
from typing import Any

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.constants import ONCOKB_API_URL, ONCOKB_LEVELS, ONCOKB_RATE_LIMIT
from variant_mcp.models.evidence import OncoKBAnnotation

logger = logging.getLogger(__name__)

ONCOKB_NO_TOKEN_MSG = (
    "OncoKB API token not configured. OncoKB provides FDA-recognized oncogenicity "
    "annotations and therapeutic levels. To enable OncoKB queries:\n"
    "1. Register for a free academic API token at https://www.oncokb.org/apiAccess\n"
    "2. Set the environment variable: ONCOKB_API_TOKEN=<your_token>\n"
    "Results from other sources (CIViC, ClinVar) are still available."
)


class OncoKBClient(BaseClient):
    """Client for the OncoKB REST API."""

    def __init__(self) -> None:
        self._token = os.environ.get("ONCOKB_API_TOKEN")
        headers = {}
        if self._token:
            headers["Authorization"] = f"Bearer {self._token}"
        super().__init__(
            base_url=ONCOKB_API_URL,
            rate_limit=ONCOKB_RATE_LIMIT,
            headers=headers,
        )

    @property
    def is_available(self) -> bool:
        return bool(self._token)

    def _check_token(self) -> str | None:
        """Return error message if token not set, else None."""
        if not self._token:
            return ONCOKB_NO_TOKEN_MSG
        return None

    async def annotate_mutation(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None = None,
    ) -> OncoKBAnnotation | str:
        """Annotate a somatic mutation by protein change.

        Returns OncoKBAnnotation on success, or error string if token missing.
        """
        token_err = self._check_token()
        if token_err:
            return token_err

        params: dict[str, str] = {
            "hugoSymbol": gene,
            "alteration": variant,
        }
        if tumor_type:
            params["tumorType"] = tumor_type

        try:
            data = await self.get_json("annotate/mutations/byProteinChange", params=params)
            return self._parse_annotation(data, gene, variant)
        except ClientError as e:
            logger.error("OncoKB annotate_mutation error: %s", e)
            raise

    async def annotate_by_hgvsg(
        self,
        hgvsg: str,
        tumor_type: str | None = None,
        ref_genome: str = "GRCh37",
    ) -> OncoKBAnnotation | str:
        """Annotate a mutation by HGVSg notation."""
        token_err = self._check_token()
        if token_err:
            return token_err

        params: dict[str, str] = {
            "hgvsg": hgvsg,
            "referenceGenome": ref_genome,
        }
        if tumor_type:
            params["tumorType"] = tumor_type

        data = await self.get_json("annotate/mutations/byHGVSg", params=params)
        return self._parse_annotation(data)

    async def get_gene_info(self, gene: str) -> dict[str, Any] | str:
        """Get gene summary and oncogene/TSG status."""
        token_err = self._check_token()
        if token_err:
            return token_err

        data = await self.get_json("genes", params={"hugoSymbol": gene})
        if isinstance(data, list) and data:
            data = data[0]
        return data

    async def get_cancer_gene_list(self) -> list[dict[str, Any]] | str:
        """Get all OncoKB cancer genes."""
        token_err = self._check_token()
        if token_err:
            return token_err

        data = await self.get_json("utils/cancerGeneList")
        return data if isinstance(data, list) else []

    async def get_all_actionable_variants(self) -> list[dict[str, Any]] | str:
        """Get all actionable therapeutic biomarkers."""
        token_err = self._check_token()
        if token_err:
            return token_err

        data = await self.get_json("utils/allActionableVariants")
        return data if isinstance(data, list) else []

    @staticmethod
    def _parse_annotation(
        data: dict[str, Any],
        gene: str | None = None,
        variant: str | None = None,
    ) -> OncoKBAnnotation:
        """Parse OncoKB annotation response."""
        mutation_effect = data.get("mutationEffect", {}) or {}
        treatments = []

        for treatment in data.get("treatments", []):
            level = treatment.get("level", "")
            drugs = [d.get("drugName", "") for d in treatment.get("drugs", [])]
            indications = []
            for indication in (
                treatment.get("levelAssociatedCancerType", {}).get("subtype", ""),
                treatment.get("levelAssociatedCancerType", {}).get("mainType", {}).get("name", ""),
            ):
                if indication:
                    indications.append(indication)

            treatments.append(
                {
                    "level": level,
                    "level_description": ONCOKB_LEVELS.get(level, ""),
                    "drugs": drugs,
                    "indications": indications,
                }
            )

        oncokb_url = None
        if gene and variant:
            oncokb_url = f"https://www.oncokb.org/gene/{gene}/{variant}"

        return OncoKBAnnotation(
            oncogenic=data.get("oncogenic"),
            mutation_effect=mutation_effect.get("description"),
            known_effect=mutation_effect.get("knownEffect"),
            highest_sensitive_level=data.get("highestSensitiveLevel"),
            highest_resistance_level=data.get("highestResistanceLevel"),
            highest_diagnostic_level=data.get("highestDiagnosticImplicationLevel"),
            highest_prognostic_level=data.get("highestPrognosticImplicationLevel"),
            treatments=treatments,
            gene_summary=data.get("geneSummary"),
            variant_summary=data.get("variantSummary"),
            oncokb_url=oncokb_url,
        )
