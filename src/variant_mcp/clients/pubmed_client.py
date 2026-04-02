"""NCBI E-utilities client for PubMed literature search."""

from __future__ import annotations

import logging
import os
from xml.etree import ElementTree

from variant_mcp.clients.base_client import BaseClient, ClientError
from variant_mcp.constants import (
    PUBMED_EFETCH_URL,
    PUBMED_ESEARCH_URL,
    PUBMED_RATE_LIMIT,
)
from variant_mcp.models.evidence import Publication, PubMedSearchResult

logger = logging.getLogger(__name__)


class PubMedClient(BaseClient):
    """Client for NCBI E-utilities (PubMed) for literature search."""

    def __init__(self) -> None:
        self._api_key = os.environ.get("NCBI_API_KEY")
        rate = 10 if self._api_key else PUBMED_RATE_LIMIT
        super().__init__(base_url="https://eutils.ncbi.nlm.nih.gov", rate_limit=rate)

    def _base_params(self) -> dict[str, str]:
        params: dict[str, str] = {}
        if self._api_key:
            params["api_key"] = self._api_key
        return params

    async def search_publications(
        self,
        gene: str,
        variant: str | None = None,
        disease: str | None = None,
        limit: int = 10,
    ) -> PubMedSearchResult:
        """Search PubMed for gene/variant/disease co-occurrence publications.

        Builds a query like: BRAF[gene] AND V600E AND melanoma
        Returns total count and list of publication summaries.

        Args:
            gene: Gene symbol, e.g. "BRAF".
            variant: Protein change, e.g. "V600E".
            disease: Disease name, e.g. "melanoma".
            limit: Maximum number of publications to return.

        Returns:
            PubMedSearchResult with total count and publications.
        """
        query_parts = [f"{gene}[gene]"]
        if variant:
            query_parts.append(variant)
        if disease:
            query_parts.append(disease)
        query = " AND ".join(query_parts)

        pmids, total_count = await self._esearch(query, retmax=limit)
        publications: list[Publication] = []
        if pmids:
            publications = await self._efetch_publications(pmids)

        return PubMedSearchResult(
            query=query,
            total_count=total_count,
            publications=publications,
        )

    async def get_publication(self, pmid: str) -> Publication:
        """Fetch publication details by PMID using efetch.

        Args:
            pmid: PubMed ID, e.g. "28138153".

        Returns:
            Publication with title, authors, abstract, etc.

        Raises:
            ClientError: If the PMID is not found or parsing fails.
        """
        publications = await self._efetch_publications([pmid])
        if not publications:
            raise ClientError(f"Publication PMID {pmid} not found.")
        return publications[0]

    async def _esearch(self, term: str, retmax: int = 10) -> tuple[list[str], int]:
        """Run esearch on PubMed and return (list of PMIDs, total count)."""
        params = self._base_params()
        params.update(
            {
                "db": "pubmed",
                "retmode": "json",
                "term": term,
                "retmax": str(retmax),
                "sort": "relevance",
            }
        )
        response = await self._request("GET", PUBMED_ESEARCH_URL, params=params)
        data = response.json()
        result = data.get("esearchresult", {})

        if "ERROR" in result:
            raise ClientError(f"PubMed esearch error: {result['ERROR']}")

        total_count = int(result.get("count", 0))
        pmids: list[str] = result.get("idlist", [])
        return pmids, total_count

    async def _efetch_publications(self, pmids: list[str]) -> list[Publication]:
        """Fetch publication details via efetch XML and parse results."""
        params = self._base_params()
        params.update(
            {
                "db": "pubmed",
                "id": ",".join(pmids),
                "rettype": "abstract",
                "retmode": "xml",
            }
        )
        response = await self._request("GET", PUBMED_EFETCH_URL, params=params)
        return self._parse_efetch_xml(response.text)

    @staticmethod
    def _parse_efetch_xml(xml_text: str) -> list[Publication]:
        """Parse PubMed efetch XML into a list of Publication models."""
        try:
            root = ElementTree.fromstring(xml_text)
        except ElementTree.ParseError as e:
            raise ClientError(f"Failed to parse PubMed XML: {e}") from e

        publications: list[Publication] = []
        for article_el in root.findall(".//PubmedArticle"):
            pub = _parse_pubmed_article(article_el)
            if pub is not None:
                publications.append(pub)

        return publications


def _parse_pubmed_article(article_el: ElementTree.Element) -> Publication | None:
    """Parse a single PubmedArticle XML element into a Publication."""
    medline = article_el.find("MedlineCitation")
    if medline is None:
        return None

    # PMID
    pmid_el = medline.find("PMID")
    if pmid_el is None or not pmid_el.text:
        return None
    pmid = pmid_el.text

    article = medline.find("Article")
    if article is None:
        return Publication(pmid=pmid, pubmed_url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")

    # Title
    title_el = article.find("ArticleTitle")
    title = title_el.text if title_el is not None and title_el.text else None

    # Authors
    authors: list[str] = []
    author_list = article.find("AuthorList")
    if author_list is not None:
        for author_el in author_list.findall("Author"):
            last = author_el.find("LastName")
            initials = author_el.find("Initials")
            if last is not None and last.text:
                name = last.text
                if initials is not None and initials.text:
                    name = f"{name} {initials.text}"
                authors.append(name)

    # Journal and year
    journal_el = article.find("Journal")
    journal: str | None = None
    year: str | None = None
    if journal_el is not None:
        journal_title = journal_el.find("Title")
        if journal_title is not None and journal_title.text:
            journal = journal_title.text
        pub_date = journal_el.find("JournalIssue/PubDate")
        if pub_date is not None:
            year_el = pub_date.find("Year")
            if year_el is not None and year_el.text:
                year = year_el.text
            elif pub_date.find("MedlineDate") is not None:
                md = pub_date.find("MedlineDate")
                if md is not None and md.text:
                    year = md.text[:4]

    # Abstract
    abstract: str | None = None
    abstract_el = article.find("Abstract")
    if abstract_el is not None:
        parts: list[str] = []
        for text_el in abstract_el.findall("AbstractText"):
            label = text_el.get("Label")
            text = "".join(text_el.itertext())
            if text:
                if label:
                    parts.append(f"{label}: {text}")
                else:
                    parts.append(text)
        if parts:
            abstract = " ".join(parts)

    # MeSH terms
    mesh_terms: list[str] = []
    mesh_list = medline.find("MeshHeadingList")
    if mesh_list is not None:
        for mesh_heading in mesh_list.findall("MeshHeading"):
            descriptor = mesh_heading.find("DescriptorName")
            if descriptor is not None and descriptor.text:
                mesh_terms.append(descriptor.text)

    # DOI
    doi: str | None = None
    article_id_list = article_el.find("PubmedData/ArticleIdList")
    if article_id_list is not None:
        for aid in article_id_list.findall("ArticleId"):
            if aid.get("IdType") == "doi" and aid.text:
                doi = aid.text
                break

    return Publication(
        pmid=pmid,
        title=title,
        authors=authors,
        journal=journal,
        year=year,
        abstract=abstract,
        mesh_terms=mesh_terms,
        doi=doi,
        pubmed_url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    )
