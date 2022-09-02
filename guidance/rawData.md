# Raw sequence data

Dear GLOMICON partners,

Below, please find some supporting guidance on how to provide your raw sequence data outputs in preparation for the cross-observatory analysis.

---

## Data Format
Please use the standard data formats for raw sequencing data: FASTA or FASTQ.

## Hosting
All files will need permanent identifiers (PIDs). Ideally, the data should be uploaded to one of the INSDC resources (ENA, GenBank, DDBJ), however, until the publishing of this project, the data can also be hosted on restricted-access servers. Possible options would be to host the data on zenodo and restrict access, or to use your institutional server with long term archiving.

## Sharing
Once you have found a long term and stable home for your raw data, please share the PIDs to the data in 
1) JSON-LD following the format as given in [LINK], with one JSON-LD per multi-fastq file
2) a csv file which we will then put into the JSON-LD format

Please note that a minimal metadata and provenance information should accompany the links to every raw sequence file to allow the clear identification of the linked out dataset.
This metadata is included in the JSON-LD file using schema.org, Dublin Core, or Darwin Core semantics for keys and consists of the following:
- name
- description
- url
- keywords
- provider
- measurementTechnique
- variableMeasured
- temporalCoverage
- spatialCoverage
- license
- version
- includedInDataCatalog (to be filled once the data is openly available on the INSDC)
