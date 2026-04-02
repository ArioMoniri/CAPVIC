"""MyVariant.info client for aggregated in-silico prediction scores."""

from __future__ import annotations

import logging
import sys
from typing import Any

from variant_mcp.clients.base_client import BaseClient
from variant_mcp.constants import MYVARIANT_API_URL, MYVARIANT_RATE_LIMIT
from variant_mcp.models.evidence import InSilicoPredictions

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stderr))


class MyVariantClient(BaseClient):
    """Client for MyVariant.info API for aggregated in-silico predictions."""

    _PREDICTION_FIELDS = (
        "dbnsfp.sift,dbnsfp.polyphen2,dbnsfp.revel,"
        "dbnsfp.cadd,dbnsfp.gerp,dbnsfp.alphamissense,dbnsfp.phylop"
    )

    def __init__(self) -> None:
        super().__init__(
            base_url=MYVARIANT_API_URL,
            rate_limit=MYVARIANT_RATE_LIMIT,
        )

    async def get_predictions(self, hgvs_id: str) -> InSilicoPredictions:
        """Get aggregated prediction scores for a variant.

        Args:
            hgvs_id: HGVS identifier, e.g. ``chr7:g.140453136A>T``.
        """
        params = {"fields": self._PREDICTION_FIELDS}
        data = await self.get_json(f"variant/{hgvs_id}", params=params)
        return self._parse_predictions(data)

    async def search_variant(self, gene: str, variant: str) -> InSilicoPredictions | None:
        """Search by gene + variant protein change.

        Args:
            gene: Gene symbol, e.g. ``BRAF``.
            variant: Protein change, e.g. ``V600E``.
        """
        params = {
            "q": f"dbnsfp.genename:{gene} AND dbnsfp.aa.change:{variant}",
            "fields": self._PREDICTION_FIELDS,
            "size": "1",
        }
        data = await self.get_json("query", params=params)
        hits: list[dict[str, Any]] = data.get("hits", [])
        if not hits:
            return None
        return self._parse_predictions(hits[0])

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_float(value: Any) -> float | None:
        """Safely extract a float from potentially nested/listed dbNSFP data."""
        if value is None:
            return None
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, list):
            for v in value:
                try:
                    return float(v)
                except (TypeError, ValueError):
                    continue
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _safe_str(value: Any) -> str | None:
        """Safely extract a string from potentially nested/listed dbNSFP data."""
        if value is None:
            return None
        if isinstance(value, str):
            return value
        if isinstance(value, list):
            for v in value:
                if isinstance(v, str) and v not in (".", ""):
                    return v
            return None
        return str(value)

    def _parse_predictions(self, data: dict[str, Any]) -> InSilicoPredictions:
        """Parse dbNSFP prediction data into an InSilicoPredictions model."""
        dbnsfp: dict[str, Any] = data.get("dbnsfp", {})
        if not dbnsfp:
            return InSilicoPredictions()

        sift = dbnsfp.get("sift", {})
        polyphen2 = dbnsfp.get("polyphen2", {})
        cadd = dbnsfp.get("cadd", {})
        gerp = dbnsfp.get("gerp++", dbnsfp.get("gerp", {}))
        am = dbnsfp.get("alphamissense", {})
        phylop = dbnsfp.get("phylop", {})
        revel = dbnsfp.get("revel", {})

        # Handle both dict and direct-value formats
        sift_score = self._safe_float(sift.get("score") if isinstance(sift, dict) else sift)
        sift_pred = self._safe_str(sift.get("pred") if isinstance(sift, dict) else None)
        pp2_score = self._safe_float(
            polyphen2.get("hdiv", {}).get("score") if isinstance(polyphen2, dict) else polyphen2
        )
        pp2_pred = self._safe_str(
            polyphen2.get("hdiv", {}).get("pred") if isinstance(polyphen2, dict) else None
        )
        revel_score = self._safe_float(revel.get("score") if isinstance(revel, dict) else revel)
        cadd_phred = self._safe_float(cadd.get("phred") if isinstance(cadd, dict) else cadd)
        am_score = self._safe_float(am.get("score") if isinstance(am, dict) else am)
        am_pred = self._safe_str(am.get("pred") if isinstance(am, dict) else None)
        gerp_score = self._safe_float(gerp.get("rs") if isinstance(gerp, dict) else gerp)
        phylop_score = self._safe_float(
            phylop.get("100way_vertebrate", {}).get("score") if isinstance(phylop, dict) else phylop
        )

        # Compute consensus
        damaging = 0
        benign = 0
        total = 0

        if sift_score is not None:
            total += 1
            if sift_score < 0.05:
                damaging += 1
            else:
                benign += 1

        if pp2_score is not None:
            total += 1
            if pp2_score > 0.85:
                damaging += 1
            elif pp2_score < 0.15:
                benign += 1

        if revel_score is not None:
            total += 1
            if revel_score > 0.5:
                damaging += 1
            else:
                benign += 1

        if cadd_phred is not None:
            total += 1
            if cadd_phred > 20:
                damaging += 1
            else:
                benign += 1

        if am_score is not None:
            total += 1
            if am_pred and "pathogenic" in am_pred.lower():
                damaging += 1
            elif am_pred and "benign" in am_pred.lower():
                benign += 1

        if gerp_score is not None:
            total += 1
            if gerp_score > 2:
                damaging += 1
            else:
                benign += 1

        consensus: str | None = None
        if total > 0:
            if damaging > benign:
                consensus = "Damaging"
            elif benign > damaging:
                consensus = "Benign"
            else:
                consensus = "Mixed"

        return InSilicoPredictions(
            sift_score=sift_score,
            sift_prediction=sift_pred,
            polyphen2_score=pp2_score,
            polyphen2_prediction=pp2_pred,
            revel_score=revel_score,
            cadd_phred=cadd_phred,
            alphamissense_score=am_score,
            alphamissense_prediction=am_pred,
            gerp_score=gerp_score,
            phylop_score=phylop_score,
            consensus=consensus,
            damaging_count=damaging,
            benign_count=benign,
            total_predictors=total,
        )
