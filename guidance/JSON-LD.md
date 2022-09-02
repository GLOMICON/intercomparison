🚧 **Please note that this page is currently under construction and may be subject to change** 🚧

Dear GLOMICON partners,

Below, please find some supporting guidance on how to provide JSON-LD files linking out to all the other resources (raw data, ASV tables, metadata, protocol). Please note that this step is optional for the partners, and will otherwise be provided through the AWI partner.

---

## Formatting
We will be following the guidance from the [Ocean InfoHub Project](https://book.oceaninfohub.org/index.html). Please have a look at their introduction to [Structured Data on the Web](https://book.oceaninfohub.org/content.html#) for some background information.

As described, we will use the JSON-LD format in combination with the [schema.org](http://schema.org/) vocabulary. 
The Ocean InfoHub project provides thematic patterns for [Documents](https://book.oceaninfohub.org/thematics/docs/README.html#) (for protocols) and for [Datasets](https://book.oceaninfohub.org/thematics/dataset/index.html) (for metadata, raw data, ASV tables), which work well for our purposes, and which we will thus be using.

The schema.org vocabularies that will be most essential are
- [schema.org Dataset markup](https://schema.org/Dataset) for raw data, metadata, ASV tables
- [schema.org CreativeWork markup](https://schema.org/CreativeWork) for protocols

Example JSON-LD files are provided for each collection of (meta)data in the respective folders in this repository. These can be copied and adapted to hold the inforamtion for all partners.

The [JSON-LD Playground](https://json-ld.org/playground/) is a great tool to check if the syntax is correct and for visualisation.



## Hosting
Please follow the following naming convention for the JSON-LD files: 
- 

The JSON-LD files will be hosted in the Github repository:
- folder name: will host all JSON-LD files including information on the metadata
- folder name: will host all JSON-LD files including information on the protocols 
- folder name: will host all JSON-LD files including information on the raw data 
- folder name: will host all JSON-LD files including information on the ASV tables

## Sharing
As the JSON-LD files will be hosted in this Github repository, no additional steps will need to be taken to share them with the partners.
