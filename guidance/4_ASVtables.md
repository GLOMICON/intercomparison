Dear GLOMICON partners,

Below, please find some supporting guidance on how to provide your raw sequence data outputs in preparation for the cross-observatory analysis.

Please let us know if you need any support. If needed, we would also be happy to offer a session where we go through an example of sharing the ASV tables via BIOM files.

---

## Data format
Please use the GSC-endorsed file format for sharing ASV tables: the BIOM format. This format allows the easy exchange and space-efficient storage of large sparse tables.
Note, that this data format can additionally include metadata tables and taxonomic information. If possible, please provide all three in this

Please see the [BIOM website](http://biom-format.org/index.html) for detailed guidance on how to use the file format, and the accompanying [publication](https://doi.org/10.1186%2F2047-217X-1-7) for more background information. 

## Hosting
All files will need permanent identifiers (PIDs). Ideally, the data should be uploaded to a dedicated stable and trusted repository, however, no such repository exists yet. Thus, please host the data on zenodo (restrict access as needed), or use your institutional server with long term archiving.

## Sharing
Once you have found a long term and stable home for your ASV tables, please share the PIDs to them in a csv file either in the dedicated folder in this repository or via email.

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

Please see [here](https://github.com/GLOMICON/intercomparison/blob/main/ASVtables/ASVtableExample.json) for an example of this metadata
