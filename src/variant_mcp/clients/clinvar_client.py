"""NCBI E-utilities client for ClinVar."""

from __future__ import annotations

import logging
import os
from typing import Any
from xml.etree import ElementTree

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.constants import (
    CLINVAR_BASE_URL,
    CLINVAR_EFETCH_URL,
    CLINVAR_ESEARCH_URL,
    CLINVAR_ESUMMARY_URL,
    CLINVAR_RATE_LIMIT,
    CLINVAR_REVIEW_STARS,
)
from variant_mcp.models.evidence import ClinVarVariant

logger = logging.getLogger(__name__)


class ClinVarClient(BaseClient):
    """Client for NCBI E-utilities (ClinVar)."""

    def __init__(self) -> None:
        self._api_key = os.environ.get("NCBI_API_KEY")
        rate = 10 if self._api_key else CLINVAR_RATE_LIMIT
        super().__init__(base_url="https://eutils.ncbi.nlm.nih.gov", rate_limit=rate)

    def _base_params(self) -> dict[str, str]:
        params: dict[str, str] = {"db": "clinvar", "retmode": "json"}
        if self._api_key:
            params["api_key"] = self._api_key
        return params

    async def search_variants(
        self,
        *,
        gene: str | None = None,
        variant_name: str | None = None,
        clinical_significance: str | None = None,
        review_status: str | None = None,
        disease: str | None = None,
        limit: int = 25,
    ) -> list[ClinVarVariant]:
        """Search ClinVar and return parsed variant summaries."""
        query_parts: list[str] = []
        if gene:
            query_parts.append(f"{gene}[gene]")
        if variant_name:
            query_parts.append(f"{variant_name}[variant name]")
        if clinical_significance:
            query_parts.append(f"{clinical_significance}[CLNSIG]")
        if review_status:
            query_parts.append(f"{review_status}[RVSTAT]")
        if disease:
            query_parts.append(f"{disease}[Disease/Phenotype]")

        if not query_parts:
            raise ClientError(
                "At least one search parameter is required (gene, variant_name, "
                "clinical_significance, disease)."
            )

        term = " AND ".join(query_parts)
        ids = await self._esearch(term, retmax=limit)
        if not ids:
            return []
        return await self._esummary_variants(ids)

    async def get_variant_detail(self, variation_id: int) -> ClinVarVariant:
        """Get full ClinVar variant record by variation ID using efetch XML."""
        params = {
            "db": "clinvar",
            "id": str(variation_id),
            "rettype": "variation",
        }
        if self._api_key:
            params["api_key"] = self._api_key

        response = await self._request("GET", CLINVAR_EFETCH_URL, params=params)
        return self._parse_efetch_xml(response.text, variation_id)

    async def get_variant_by_rsid(self, rsid: str) -> list[ClinVarVariant]:
        """Search ClinVar by dbSNP rsID."""
        clean_rsid = rsid.lower().replace("rs", "")
        ids = await self._esearch(f"rs{clean_rsid}[dbsnp_id]")
        if not ids:
            return []
        return await self._esummary_variants(ids)

    async def get_variant_by_hgvs(self, hgvs_expression: str) -> list[ClinVarVariant]:
        """Search ClinVar by HGVS expression."""
        ids = await self._esearch(f"{hgvs_expression}[variant name]")
        if not ids:
            return []
        return await self._esummary_variants(ids)

    async def _esearch(self, term: str, retmax: int = 25) -> list[str]:
        """Run esearch and return list of IDs."""
        params = self._base_params()
        params.update({"term": term, "retmax": str(retmax)})
        response = await self._request("GET", CLINVAR_ESEARCH_URL, params=params)
        data = response.json()
        result = data.get("esearchresult", {})
        if "ERROR" in result:
            raise ClientError(f"ClinVar esearch error: {result['ERROR']}")
        return result.get("idlist", [])  # type: ignore[no-any-return]

    async def _esummary_variants(self, ids: list[str]) -> list[ClinVarVariant]:
        """Run esummary for given IDs and parse to ClinVarVariant list."""
        params = self._base_params()
        params["id"] = ",".join(ids)
        response = await self._request("GET", CLINVAR_ESUMMARY_URL, params=params)
        data = response.json()
        result = data.get("result", {})
        uid_list = result.get("uids", [])

        variants = []
        for uid in uid_list:
            doc = result.get(uid, {})
            if not doc:
                continue
            variants.append(self._parse_esummary_doc(doc, uid))
        return variants

    @staticmethod
    def _parse_esummary_doc(doc: dict[str, Any], uid: str) -> ClinVarVariant:
        """Parse an esummary document into a ClinVarVariant."""
        genes_list = doc.get("genes", [])
        gene_names = [g.get("symbol", "") for g in genes_list if g.get("symbol")]

        clin_sig = doc.get("clinical_significance", {})
        if isinstance(clin_sig, dict):
            significance = clin_sig.get("description", "")
            review_status = clin_sig.get("review_status", "")
            last_eval = clin_sig.get("last_evaluated", "")
        else:
            significance = str(clin_sig)
            review_status = ""
            last_eval = ""

        review_stars = CLINVAR_REVIEW_STARS.get(review_status.lower(), "")

        variation_id = int(uid)
        return ClinVarVariant(
            variation_id=variation_id,
            title=doc.get("title", ""),
            clinical_significance=significance,
            review_status=review_status,
            review_stars=review_stars,
            genes=gene_names,
            variation_type=doc.get("variation_type", ""),
            protein_change=doc.get("protein_change", ""),
            conditions=[
                t.get("trait_name", "") for t in doc.get("trait_set", []) if t.get("trait_name")
            ],
            last_evaluated=last_eval,
            clinvar_url=f"{CLINVAR_BASE_URL}/{variation_id}",
        )

    @staticmethod
    def _parse_efetch_xml(xml_text: str, variation_id: int) -> ClinVarVariant:
        """Parse ClinVar efetch XML response."""
        try:
            root = ElementTree.fromstring(xml_text)
        except ElementTree.ParseError as e:
            raise ClientError(f"Failed to parse ClinVar XML: {e}") from e

        variant = ClinVarVariant(
            variation_id=variation_id,
            clinvar_url=f"{CLINVAR_BASE_URL}/{variation_id}",
        )

        # Parse VariationArchive
        va = root.find(".//VariationArchive")
        if va is not None:
            variant.title = va.get("VariationName", "")
            variant.variation_type = va.get("VariationType", "")

        # Parse InterpretedRecord
        ir = root.find(".//InterpretedRecord")
        if ir is None:
            ir = root.find(".//ClassifiedRecord")

        if ir is not None:
            # Clinical significance
            interp = ir.find(".//Interpretation") or ir.find(".//Classification")
            if interp is not None:
                desc = interp.find("Description")
                if desc is not None and desc.text:
                    variant.clinical_significance = desc.text
                review = interp.find("ReviewStatus")
                if review is not None and review.text:
                    variant.review_status = review.text
                    variant.review_stars = CLINVAR_REVIEW_STARS.get(review.text.lower(), "")
                last_eval_el = interp.find("DateLastEvaluated")
                if last_eval_el is not None and last_eval_el.text:
                    variant.last_evaluated = last_eval_el.text

            # HGVS expressions
            hgvs_list = []
            for hgvs_el in ir.findall(".//HGVSlist/HGVS") + ir.findall(".//HGVSExpression"):
                text = hgvs_el.text or hgvs_el.get("Expression", "")
                if text:
                    hgvs_list.append(text)
            variant.hgvs_expressions = hgvs_list

            # Genes
            gene_names = []
            for gene_el in ir.findall(".//GeneList/Gene"):
                symbol = gene_el.get("Symbol")
                if symbol:
                    gene_names.append(symbol)
            variant.genes = gene_names

            # Conditions
            conditions = []
            for trait in ir.findall(".//TraitSet/Trait"):
                name_el = trait.find("Name/ElementValue[@Type='Preferred']")
                if name_el is not None and name_el.text:
                    conditions.append(name_el.text)
            variant.conditions = conditions

            # Submitter classifications
            submitter_list = []
            for scv in ir.findall(".//ClinicalAssertionList/ClinicalAssertion"):
                submitter_el = scv.find(".//ClinVarAccession")
                interp_el = scv.find(".//Interpretation/Description")
                if interp_el is None:
                    interp_el = scv.find(".//Classification/Description")
                if submitter_el is not None:
                    entry = {
                        "submitter": submitter_el.get("SubmitterName", "Unknown"),
                        "accession": submitter_el.get("Accession", ""),
                        "classification": interp_el.text if interp_el is not None else "N/A",
                    }
                    submitter_list.append(entry)
            variant.submitter_classifications = submitter_list
            variant.submitter_count = len(submitter_list)

            # Check for conflicting interpretations
            classifications = {
                s["classification"] for s in submitter_list if s["classification"] != "N/A"
            }
            variant.conflicting = len(classifications) > 1

        # rsID
        for xref in root.findall(".//XRef"):
            if xref.get("DB") == "dbSNP":
                variant.rsid = f"rs{xref.get('ID', '')}"
                break

        return variant
