# Raw sequence data

Dear GLOMICON partners,

Below, please find some supporting guidance on how to provide your raw sequence data outputs in preparation for the cross-observatory analysis.

---

## Data Format
Please use the standard data formats for raw sequencing data: FASTA or FASTQ.

## Hosting
All files will need permanent identifiers (PIDs). Ideally, the data should be uploaded to one of the INSDC resources (ENA, GenBank, DDBJ), however, until the publishing of this project, the data can also be hosted on restricted-access servers. Possible options would be to host the data on zenodo (restrict access as needed), or to use your institutional server with long term archiving.

## Sharing
Once the raw data are safely stored in a long term and stable repository, please share the PIDs to the data in a csv file either in the dedicated folder in this repository or via email. 

We will then prepare JSON-LD files with the PIDs, provenance information and some minimal metadata. 

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

Please see [here](https://github.com/GLOMICON/intercomparison/blob/main/rawData/rawDataExample.json) for an example of the metadata needed.
