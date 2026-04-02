"""CIViC V2 GraphQL query strings."""


class CIViCQueries:
    """GraphQL queries for the CIViC V2 API."""

    EVIDENCE_FIELDS = """
        id
        name
        status
        evidenceType
        evidenceLevel
        evidenceDirection
        evidenceRating
        significance
        description
        therapyInteractionType
        molecularProfile {
            id
            name
            variants {
                name
                id
            }
        }
        disease {
            id
            name
            doid
        }
        therapies {
            id
            name
            ncitId
        }
        phenotypes {
            id
            name
            hpoId
        }
        source {
            id
            citation
            sourceUrl
            citationId
            sourceType
        }
    """

    SEARCH_EVIDENCE = (
        """
    query SearchEvidence(
        $diseaseName: String,
        $molecularProfileName: String,
        $therapyName: String,
        $evidenceType: EvidenceType,
        $significance: EvidenceSignificance,
        $first: Int,
        $after: String
    ) {
        evidenceItems(
            diseaseName: $diseaseName,
            molecularProfileName: $molecularProfileName,
            therapyName: $therapyName,
            evidenceType: $evidenceType,
            significance: $significance,
            status: ACCEPTED,
            first: $first,
            after: $after
        ) {
            totalCount
            pageInfo {
                hasNextPage
                endCursor
            }
            nodes {
                """
        + EVIDENCE_FIELDS
        + """
            }
        }
    }
    """
    )

    GET_GENE = """
    query GetGene($name: String!) {
        genes(name: $name) {
            nodes {
                id
                name
                officialName
                description
                entrezId
                variants {
                    totalCount
                    nodes {
                        id
                        name
                        singleVariantMolecularProfileId
                    }
                }
            }
        }
    }
    """

    GET_VARIANT = """
    query GetVariant($id: Int!) {
        variant(id: $id) {
            id
            name
            variantTypes {
                id
                name
                soid
            }
            singleVariantMolecularProfile {
                id
                name
                molecularProfileScore
            }
            gene {
                id
                name
            }
        }
    }
    """

    GET_EVIDENCE_ITEM = (
        """
    query GetEvidenceItem($id: Int!) {
        evidenceItem(id: $id) {
            """
        + EVIDENCE_FIELDS
        + """
        }
    }
    """
    )

    SEARCH_ASSERTIONS = """
    query SearchAssertions(
        $diseaseName: String,
        $molecularProfileName: String,
        $therapyName: String,
        $significance: String,
        $first: Int,
        $after: String
    ) {
        assertions(
            diseaseName: $diseaseName,
            molecularProfileName: $molecularProfileName,
            therapyName: $therapyName,
            significance: $significance,
            status: ACCEPTED,
            first: $first,
            after: $after
        ) {
            totalCount
            pageInfo {
                hasNextPage
                endCursor
            }
            nodes {
                id
                name
                assertionType
                assertionDirection
                significance
                status
                summary
                description
                nccnGuideline {
                    name
                }
                acmgCodes {
                    id
                    code
                    description
                }
                ampLevel
                molecularProfile {
                    id
                    name
                    variants {
                        name
                        id
                    }
                }
                disease {
                    id
                    name
                    doid
                }
                therapies {
                    id
                    name
                    ncitId
                }
                evidenceItems {
                    id
                    evidenceType
                    evidenceLevel
                    significance
                }
            }
        }
    }
    """

    TYPEAHEAD_GENES = """
    query TypeaheadGenes($queryTerm: String!) {
        geneTypeahead(queryTerm: $queryTerm) {
            id
            name
            entrezId
        }
    }
    """

    TYPEAHEAD_DISEASES = """
    query TypeaheadDiseases($queryTerm: String!) {
        diseaseTypeahead(queryTerm: $queryTerm) {
            id
            name
            doid
        }
    }
    """

    TYPEAHEAD_THERAPIES = """
    query TypeaheadTherapies($queryTerm: String!) {
        therapyTypeahead(queryTerm: $queryTerm) {
            id
            name
            ncitId
        }
    }
    """
