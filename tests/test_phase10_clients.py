"""Mock-based HTTP tests for Phase 10 clients (gnomAD, PubMed, UniProt, MyVariant)."""

from __future__ import annotations

import httpx
import pytest
import respx

from variant_mcp.clients.base_client import ClientError
from variant_mcp.clients.gnomad_client import GnomADClient
from variant_mcp.clients.myvariant_client import MyVariantClient
from variant_mcp.clients.pubmed_client import PubMedClient
from variant_mcp.clients.uniprot_client import UniProtClient

# ============================================================================
# gnomAD Client
# ============================================================================


class TestGnomADClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_get_variant_frequency(self):
        client = GnomADClient()
        respx.post("https://gnomad.broadinstitute.org/api").mock(
            return_value=httpx.Response(
                200,
                json={
                    "data": {
                        "variant": {
                            "variant_id": "7-140453136-A-T",
                            "rsids": ["rs113488022"],
                            "flags": [],
                            "genome": {
                                "ac": 3,
                                "an": 152312,
                                "af": 1.97e-05,
                                "ac_hom": 0,
                                "filters": [],
                                "populations": [
                                    {"id": "nfe", "ac": 2, "an": 100000, "af": 2e-05},
                                    {"id": "afr", "ac": 1, "an": 40000, "af": 2.5e-05},
                                    {"id": "eas", "ac": 0, "an": 10000, "af": 0.0},
                                ],
                            },
                            "exome": None,
                        }
                    }
                },
            )
        )
        freq = await client.get_variant_frequency("7-140453136-A-T")
        assert freq.variant_id == "7-140453136-A-T"
        assert freq.rsid == "rs113488022"
        assert freq.allele_frequency == pytest.approx(1.97e-05)
        assert freq.allele_count == 3
        assert freq.allele_number == 152312
        assert freq.homozygote_count == 0
        assert "nfe" in freq.population_frequencies
        assert freq.population_frequencies["nfe"] == pytest.approx(2e-05)
        assert freq.genome_version == "GRCh38"
        assert freq.source == "gnomAD v4"
        assert freq.filter_status == "PASS"
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_variant_not_found(self):
        client = GnomADClient()
        respx.post("https://gnomad.broadinstitute.org/api").mock(
            return_value=httpx.Response(200, json={"data": {"variant": None}})
        )
        with pytest.raises(ClientError, match="not found"):
            await client.get_variant_frequency("99-999999-A-T")
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_graphql_error(self):
        client = GnomADClient()
        respx.post("https://gnomad.broadinstitute.org/api").mock(
            return_value=httpx.Response(
                200, json={"errors": [{"message": "Invalid variant ID format"}]}
            )
        )
        with pytest.raises(ClientError, match="GraphQL error"):
            await client.get_variant_frequency("bad-id")
        await client.close()

    @pytest.mark.asyncio
    async def test_unsupported_genome_version(self):
        client = GnomADClient()
        with pytest.raises(ClientError, match="Unsupported genome version"):
            await client.get_variant_frequency("7-140453136-A-T", genome_version="GRCh99")
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_exome_fallback(self):
        """When genome data is missing, use exome data."""
        client = GnomADClient()
        respx.post("https://gnomad.broadinstitute.org/api").mock(
            return_value=httpx.Response(
                200,
                json={
                    "data": {
                        "variant": {
                            "variant_id": "7-140453136-A-T",
                            "rsids": [],
                            "flags": ["lcr"],
                            "genome": None,
                            "exome": {
                                "ac": 10,
                                "an": 250000,
                                "af": 4e-05,
                                "ac_hom": 0,
                                "filters": ["AC0"],
                                "populations": [
                                    {"id": "sas", "ac": 5, "an": 30000, "af": 1.67e-04},
                                ],
                            },
                        }
                    }
                },
            )
        )
        freq = await client.get_variant_frequency("7-140453136-A-T")
        assert freq.allele_frequency == pytest.approx(4e-05)
        assert freq.allele_count == 10
        assert "sas" in freq.population_frequencies
        assert "lcr" in freq.filter_status
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_grch37_dataset(self):
        client = GnomADClient()

        def check_request(request: httpx.Request) -> httpx.Response:
            import json

            body = json.loads(request.content)
            assert body["variables"]["datasetId"] == "gnomad_r2_1"
            return httpx.Response(
                200,
                json={
                    "data": {
                        "variant": {
                            "variant_id": "7-140453136-A-T",
                            "rsids": [],
                            "flags": [],
                            "genome": {
                                "ac": 1,
                                "an": 30000,
                                "af": 3.3e-05,
                                "ac_hom": 0,
                                "filters": [],
                                "populations": [],
                            },
                            "exome": None,
                        }
                    }
                },
            )

        respx.post("https://gnomad.broadinstitute.org/api").mock(side_effect=check_request)
        freq = await client.get_variant_frequency("7-140453136-A-T", genome_version="GRCh37")
        assert freq.genome_version == "GRCh37"
        assert freq.source == "gnomAD v2"
        await client.close()


# ============================================================================
# PubMed Client
# ============================================================================


class TestPubMedClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_search_publications(self):
        client = PubMedClient()
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi").mock(
            return_value=httpx.Response(
                200, json={"esearchresult": {"idlist": ["28138153"], "count": "1"}}
            )
        )
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi").mock(
            return_value=httpx.Response(
                200,
                text="""<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>28138153</PMID>
      <Article>
        <ArticleTitle>CIViC is a community knowledgebase.</ArticleTitle>
        <AuthorList>
          <Author><LastName>Griffith</LastName><Initials>M</Initials></Author>
          <Author><LastName>Spies</LastName><Initials>NC</Initials></Author>
        </AuthorList>
        <Journal>
          <Title>Nature Genetics</Title>
          <JournalIssue><PubDate><Year>2017</Year></PubDate></JournalIssue>
        </Journal>
        <Abstract>
          <AbstractText>CIViC is an expert-curated knowledgebase.</AbstractText>
        </Abstract>
      </Article>
      <MeshHeadingList>
        <MeshHeading><DescriptorName>Neoplasms</DescriptorName></MeshHeading>
      </MeshHeadingList>
    </MedlineCitation>
    <PubmedData>
      <ArticleIdList>
        <ArticleId IdType="doi">10.1038/ng.3774</ArticleId>
      </ArticleIdList>
    </PubmedData>
  </PubmedArticle>
</PubmedArticleSet>""",
            )
        )
        result = await client.search_publications("BRAF", variant="V600E")
        assert result.total_count == 1
        assert len(result.publications) == 1
        pub = result.publications[0]
        assert pub.pmid == "28138153"
        assert pub.title == "CIViC is a community knowledgebase."
        assert "Griffith M" in pub.authors
        assert pub.journal == "Nature Genetics"
        assert pub.year == "2017"
        assert pub.doi == "10.1038/ng.3774"
        assert "Neoplasms" in pub.mesh_terms
        assert pub.abstract is not None
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_get_publication_not_found(self):
        client = PubMedClient()
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi").mock(
            return_value=httpx.Response(
                200, text='<?xml version="1.0"?><PubmedArticleSet></PubmedArticleSet>'
            )
        )
        with pytest.raises(ClientError, match="not found"):
            await client.get_publication("99999999999")
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_search_no_results(self):
        client = PubMedClient()
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi").mock(
            return_value=httpx.Response(200, json={"esearchresult": {"idlist": [], "count": "0"}})
        )
        result = await client.search_publications("FAKEGENE123")
        assert result.total_count == 0
        assert result.publications == []
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_esearch_error(self):
        client = PubMedClient()
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi").mock(
            return_value=httpx.Response(200, json={"esearchresult": {"ERROR": "Invalid query"}})
        )
        with pytest.raises(ClientError, match="esearch error"):
            await client.search_publications("BRAF")
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_malformed_xml(self):
        client = PubMedClient()
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi").mock(
            return_value=httpx.Response(
                200, json={"esearchresult": {"idlist": ["123"], "count": "1"}}
            )
        )
        respx.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi").mock(
            return_value=httpx.Response(200, text="<not valid xml")
        )
        with pytest.raises(ClientError, match="parse"):
            await client.search_publications("BRAF")
        await client.close()


# ============================================================================
# UniProt Client
# ============================================================================


class TestUniProtClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_get_protein_features(self):
        client = UniProtClient()
        respx.get("https://rest.uniprot.org/uniprotkb/search").mock(
            return_value=httpx.Response(
                200,
                json={
                    "results": [
                        {
                            "primaryAccession": "P15056",
                            "sequence": {"length": 766},
                            "features": [
                                {
                                    "type": "Domain",
                                    "description": "Protein kinase",
                                    "location": {
                                        "start": {"value": 457},
                                        "end": {"value": 717},
                                    },
                                    "evidences": [{"source": {"name": "Pfam"}}],
                                },
                                {
                                    "type": "Region",
                                    "description": "Pkinase_Tyr",
                                    "location": {
                                        "start": {"value": 458},
                                        "end": {"value": 712},
                                    },
                                    "evidences": [],
                                },
                            ],
                        }
                    ]
                },
            )
        )
        features = await client.get_protein_features("BRAF")
        assert features.gene == "BRAF"
        assert features.uniprot_id == "P15056"
        assert features.protein_length == 766
        assert len(features.domains) == 2
        assert features.domains[0].name == "Protein kinase"
        assert features.domains[0].start_pos == 457
        assert features.domains[0].end_pos == 717
        assert features.domains[0].source == "Pfam"
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_gene_not_found(self):
        client = UniProtClient()
        respx.get("https://rest.uniprot.org/uniprotkb/search").mock(
            return_value=httpx.Response(200, json={"results": []})
        )
        features = await client.get_protein_features("FAKEGENE")
        assert features.gene == "FAKEGENE"
        assert features.uniprot_id is None
        assert features.domains == []
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_check_variant_in_domain(self):
        client = UniProtClient()
        respx.get("https://rest.uniprot.org/uniprotkb/search").mock(
            return_value=httpx.Response(
                200,
                json={
                    "results": [
                        {
                            "primaryAccession": "P15056",
                            "sequence": {"length": 766},
                            "features": [
                                {
                                    "type": "Domain",
                                    "description": "Protein kinase",
                                    "location": {
                                        "start": {"value": 457},
                                        "end": {"value": 717},
                                    },
                                    "evidences": [],
                                }
                            ],
                        }
                    ]
                },
            )
        )
        result = await client.check_variant_in_domain("BRAF", "V600E")
        assert result.in_domain is True
        assert result.evidence_for_om1 is True
        assert result.position == 600
        assert len(result.domains) == 1
        assert result.domains[0].name == "Protein kinase"
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_variant_outside_domain(self):
        client = UniProtClient()
        respx.get("https://rest.uniprot.org/uniprotkb/search").mock(
            return_value=httpx.Response(
                200,
                json={
                    "results": [
                        {
                            "primaryAccession": "P15056",
                            "sequence": {"length": 766},
                            "features": [
                                {
                                    "type": "Domain",
                                    "description": "Protein kinase",
                                    "location": {
                                        "start": {"value": 457},
                                        "end": {"value": 717},
                                    },
                                    "evidences": [],
                                }
                            ],
                        }
                    ]
                },
            )
        )
        result = await client.check_variant_in_domain("BRAF", "A100T")
        assert result.in_domain is False
        assert result.evidence_for_om1 is False
        assert result.position == 100
        await client.close()


# ============================================================================
# MyVariant.info Client
# ============================================================================


class TestMyVariantClient:
    @respx.mock
    @pytest.mark.asyncio
    async def test_get_predictions(self):
        client = MyVariantClient()
        respx.get("https://myvariant.info/v1/variant/chr7:g.140453136A%3ET").mock(
            return_value=httpx.Response(
                200,
                json={
                    "dbnsfp": {
                        "sift": {"score": 0.0, "pred": "D"},
                        "polyphen2": {"hdiv": {"score": 1.0, "pred": "D"}},
                        "revel": {"score": 0.95},
                        "cadd": {"phred": 33.0},
                        "alphamissense": {"score": 0.99, "pred": "likely_pathogenic"},
                        "gerp++": {"rs": 5.48},
                        "phylop": {"100way_vertebrate": {"score": 9.3}},
                    }
                },
            )
        )
        preds = await client.get_predictions("chr7:g.140453136A>T")
        assert preds.sift_score == pytest.approx(0.0)
        assert preds.sift_prediction == "D"
        assert preds.polyphen2_score == pytest.approx(1.0)
        assert preds.revel_score == pytest.approx(0.95)
        assert preds.cadd_phred == pytest.approx(33.0)
        assert preds.alphamissense_score == pytest.approx(0.99)
        assert preds.gerp_score == pytest.approx(5.48)
        assert preds.phylop_score == pytest.approx(9.3)
        assert preds.consensus == "Damaging"
        assert preds.damaging_count >= 5
        assert preds.total_predictors >= 6
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_search_variant(self):
        client = MyVariantClient()
        respx.get("https://myvariant.info/v1/query").mock(
            return_value=httpx.Response(
                200,
                json={
                    "hits": [
                        {
                            "dbnsfp": {
                                "sift": {"score": 0.5, "pred": "T"},
                                "polyphen2": {"hdiv": {"score": 0.05, "pred": "B"}},
                                "revel": {"score": 0.1},
                                "cadd": {"phred": 8.0},
                            }
                        }
                    ]
                },
            )
        )
        preds = await client.search_variant("BRAF", "V600E")
        assert preds is not None
        assert preds.consensus == "Benign"
        assert preds.benign_count >= 3
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_search_variant_not_found(self):
        client = MyVariantClient()
        respx.get("https://myvariant.info/v1/query").mock(
            return_value=httpx.Response(200, json={"hits": []})
        )
        result = await client.search_variant("FAKEGENE", "X999Y")
        assert result is None
        await client.close()

    @respx.mock
    @pytest.mark.asyncio
    async def test_empty_dbnsfp(self):
        client = MyVariantClient()
        respx.get("https://myvariant.info/v1/variant/chr1:g.1A%3ET").mock(
            return_value=httpx.Response(200, json={})
        )
        preds = await client.get_predictions("chr1:g.1A>T")
        assert preds.consensus is None
        assert preds.total_predictors == 0
        await client.close()

    def test_safe_float_edge_cases(self):
        client = MyVariantClient()
        assert client._safe_float(None) is None
        assert client._safe_float(0.5) == pytest.approx(0.5)
        assert client._safe_float(1) == pytest.approx(1.0)
        assert client._safe_float("0.3") == pytest.approx(0.3)
        assert client._safe_float([0.1, 0.2]) == pytest.approx(0.1)
        assert client._safe_float([".", "0.5"]) == pytest.approx(0.5)
        assert client._safe_float("not_a_number") is None
        assert client._safe_float([]) is None

    def test_safe_str_edge_cases(self):
        client = MyVariantClient()
        assert client._safe_str(None) is None
        assert client._safe_str("deleterious") == "deleterious"
        assert client._safe_str([".", "", "tolerated"]) == "tolerated"
        assert client._safe_str(42) == "42"
        assert client._safe_str([".", ""]) is None
