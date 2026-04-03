"""Tests for API clients using mocked HTTP responses."""

import httpx
import pytest
import respx

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.clients.civic_client import CIViCClient
from variant_mcp.clients.clinvar_client import ClinVarClient
from variant_mcp.clients.oncokb_client import OncoKBClient


class TestBaseClient:
    @pytest.fixture
    def client(self):
        return BaseClient("https://example.com/api", rate_limit=10)

    @respx.mock
    @pytest.mark.asyncio
    async def test_get_json_success(self, client):
        respx.get("https://example.com/api/test").mock(
            return_value=httpx.Response(200, json={"result": "ok"})
        )
        result = await client.get_json("test")
        assert result == {"result": "ok"}
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_get_404_raises(self, client):
        respx.get("https://example.com/api/missing").mock(
            return_value=httpx.Response(404, text="Not found")
        )
        with pytest.raises(ClientError):
            await client.get_json("missing")
        await client.close()


class TestCIViCClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_search_evidence_parsed(self):
        client = CIViCClient()
        respx.post("https://civicdb.org/api/graphql").mock(
            return_value=httpx.Response(
                200,
                json={
                    "data": {
                        "evidenceItems": {
                            "totalCount": 1,
                            "pageInfo": {"hasNextPage": False, "endCursor": None},
                            "nodes": [
                                {
                                    "id": 42,
                                    "name": "EID42",
                                    "status": "ACCEPTED",
                                    "evidenceType": "PREDICTIVE",
                                    "evidenceLevel": "A",
                                    "evidenceDirection": "SUPPORTS",
                                    "evidenceRating": 5,
                                    "significance": "SENSITIVITY",
                                    "description": "Test evidence",
                                    "therapyInteractionType": None,
                                    "molecularProfile": {
                                        "id": 1,
                                        "name": "BRAF V600E",
                                        "variants": [{"name": "V600E", "id": 12}],
                                    },
                                    "disease": {"id": 1, "name": "Melanoma", "doid": "1909"},
                                    "therapies": [
                                        {"id": 1, "name": "Vemurafenib", "ncitId": "C64768"}
                                    ],
                                    "phenotypes": [],
                                    "source": {
                                        "id": 1,
                                        "citation": "Chapman 2011",
                                        "sourceUrl": None,
                                        "citationId": "21639808",
                                        "sourceType": "PUBMED",
                                    },
                                }
                            ],
                        }
                    }
                },
            )
        )
        items = await client.search_evidence_parsed(gene="BRAF", variant="V600E")
        assert len(items) == 1
        assert items[0].id == 42
        assert items[0].evidence_level == "A"
        assert "Vemurafenib" in items[0].therapies
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_search_typeahead_gene(self):
        client = CIViCClient()
        respx.post("https://civicdb.org/api/graphql").mock(
            return_value=httpx.Response(
                200,
                json={
                    "data": {
                        "featureTypeahead": [
                            {"id": 5, "name": "BRAF"},
                            {"id": 6, "name": "BRCA1"},
                        ]
                    }
                },
            )
        )
        results = await client.search_typeahead("gene", "BRA")
        assert len(results) == 2
        assert results[0]["name"] == "BRAF"
        await client.close()


class TestClinVarClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_search_variants(self):
        client = ClinVarClient()

        # Mock esearch
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi").mock(
            return_value=httpx.Response(
                200, json={"esearchresult": {"idlist": ["13961"], "count": "1"}}
            )
        )
        # Mock esummary
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi").mock(
            return_value=httpx.Response(
                200,
                json={
                    "result": {
                        "uids": ["13961"],
                        "13961": {
                            "uid": "13961",
                            "title": "NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)",
                            "clinical_significance": {
                                "description": "Pathogenic/Likely pathogenic",
                                "review_status": "reviewed by expert panel",
                                "last_evaluated": "2024-01-15",
                            },
                            "genes": [{"symbol": "BRAF"}],
                            "variation_type": "single nucleotide variant",
                            "protein_change": "V600E",
                            "trait_set": [{"trait_name": "Melanoma"}],
                        },
                    }
                },
            )
        )

        results = await client.search_variants(gene="BRAF", variant_name="V600E")
        assert len(results) == 1
        assert results[0].variation_id == 13961
        assert "Pathogenic" in results[0].clinical_significance
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_search_no_results(self):
        client = ClinVarClient()
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi").mock(
            return_value=httpx.Response(200, json={"esearchresult": {"idlist": [], "count": "0"}})
        )
        results = await client.search_variants(gene="FAKEGENE")
        assert results == []
        await client.close()


class TestOncoKBClient:
    def test_no_token_not_available(self, monkeypatch):
        monkeypatch.delenv("ONCOKB_API_TOKEN", raising=False)
        client = OncoKBClient()
        assert not client.is_available

    @pytest.mark.asyncio
    async def test_annotate_without_token(self, monkeypatch):
        monkeypatch.delenv("ONCOKB_API_TOKEN", raising=False)
        client = OncoKBClient()
        result = await client.annotate_mutation("BRAF", "V600E")
        assert isinstance(result, str)
        assert "token" in result.lower()
        await client.close()
