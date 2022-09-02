Dear GLOMICON partners, 

Below, please find some supporting guidance to help us get our metadata in order, in preparation for the overarching analysis and submission to the INSDC. 

For this project, we’d like a coherent set of metadata to facilitate comparative work and set an example for the wider community. To accomplish this we’ll be closely following the Minimum Information about any (x) Sequence (MIxS) standard. 

We will support this process and are happy to validate the metadata records and provide feedback.

---

1. As a first step, please email the metadata as simple CSVs / tables / spreadsheets to Raïssa Meyer, these will be added to GitHub as they are. 

2. In a second step, please standardise the metadata according to the Minimum Information about any (x) Sequence (MIxS) checklists (see guidance below). Please do so either in this GitHub repository (version control included), or do so on your local machine and email a simple CSV / table to Raïssa Meyer. As noted in the [guidance on ASV tables](https://github.com/GLOMICON/intercomparison/blob/main/guidance/ASVtables.md#data-format), the standardised metadata can also be added to the BIOM file.

3. In a third step, the standardised metadata will be converted into JSON-LD by us. 

---

**For the MIxS-compliant standardisation process, please see the guidance below:**

[Understanding the MIxS checklist](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#understanding-the-mixs-checklist)
- [Tab 1: README](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#tab-1-readme)
- [Tab 2: MIxS core](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#tab-2-mixs-core)
- [Tab 3: water environmental_packages](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#tab-3-water-environmental_packages)
  
[Using the MIxS checklist](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#using-the-mixs-checklist)
- [Value syntax and units](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#value-syntax-and-units)
- [Missing values](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#missing-values)
- [Ontology terms](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#ontology-terms-relevant-only-for-mixs-core) (relevant only for MIxS core)
  - [experimental_factor](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#experimental_factor)
  - [geo_loc_name](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#geo_loc_name)
  - [env_[broad_scale/local_scale/medium]](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#env_broad_scalelocal_scalemedium)
  - [samp_mat_processing](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#samp_mat_processing)

--- 

## Understanding the MIxS checklist
Our metadata guidelines are based on the MIxS checklist. To support your work, we have prepared metadata tables for each partner organisation in which relevant fields for this sample exchange are highlighted and which contain examples appropriate to your work (please see [here](https://drive.google.com/drive/folders/16e5EKksO6G4TSGs-XX-KoWC6z5QGGHiO?usp=sharing) and use the file which has your institution’s name in the title). 

The metadata spreadsheet contains three tabs:
1. [Tab 1: README](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#tab-1-readme)
2. [Tab 2: MIxS core](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#tab-2-mixs-core)
3. [Tab 3: water environmental_packages](https://github.com/GLOMICON/intercomparison/blob/main/guidance/metadata.md#tab-3-water-environmental_packages)

### Tab 1: README
The README provides general information on how the spreadsheets are organised, as well as some general guidance to approach and use the standard.

### Tab 2: MIxS core
This tab contains a table with names, definitions, syntax and unit specifications, and additional guidance on the core MIxS fields. 

### Tab 3: water environmental_packages
This tab contains a table with names, definitions and additional guidance from the MIxS environmental package **water**, which we will use to further describe the environment we sampled. 

--- 

In Tab 2 and 3, fields are highlighted in the following manner:
- Fields marked in blue are relevant for **Phase 1** and are mandatory or highly recommended based on our MIxS extension and our sample exchange
  - Phase 1 spans the sample collection up to sequence analysis 
- Fields marked in green are relevant for **Phase 2** and are mandatory or highly recommended based on our MIxS extension and our sample exchange
  - We will be entering Phase 2 when we are making our sequence data, protocols, etc. publicly available

Please also consider if any of the unmarked fields may be relevant to you. This may be especially important for the environmental package **_water_** as the amount of data recorded may vary between partners.

---

## Using the MIxS checklist

We’ve structured the spreadsheet for MIxS core and the MIxS environmental_package to have mandatory MIxS fields at the top, followed by the remaining fields. 

When entering entries in the sheets, please consider the following:

### Value syntax and units
**Please follow the value syntax** (given in tab [MIxS, column E](https://docs.google.com/spreadsheets/d/1iT2DBokrXKkf25EWrUSNeXJXXd1e46um/edit#gid=937998399&range=E:E)) **and the preferred unit** (if noted in [MIxS, column S](https://docs.google.com/spreadsheets/d/1iT2DBokrXKkf25EWrUSNeXJXXd1e46um/edit#gid=937998399&range=S:S))

### Missing values
**If a value is missing, please follow the [INSDC missing value vocabulary](https://ena-docs.readthedocs.io/en/latest/submit/samples/missing-values.html)**
- not applicable - information is inappropriate to report, can indicate that the standard itself fails to model or represent the information appropriately
- not collected - information of an expected format was not given because it has not been collected
- not provided - information of an expected format was not given, a value may be given at the later stage
- restricted access - information exists but can not be released openly because of privacy concerns

### Ontology terms (relevant only for MIxS core)
**Please be sure to use ontology terms wherever requested.**

Based on the information you have already shared, we have prepared lists of potentially relevant ontology terms (see below). 

Please let us know if none of these terms are applicable/you have done something else. [Here](https://github.com/EnvironmentOntology/envo/wiki/ENVO-annotations-for-MIxS-v5) is also some dedicated guidance on using ENVO terms to fill in metadata. If you don’t find the term you are looking for in ENVO, please let us know and we’ll add them.

If you would like additional help browsing ontologies please contact us. 

#### experimental_factor

Consider using the following [OBI](https://www.ebi.ac.uk/ols/ontologies/obi) or [EFO](https://www.ebi.ac.uk/ols/ontologies/efo) terms: 
- [methodological variation design](https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0004669&viewMode=All&siblings=true) [EFO:EFO_0004669]
  - [optimization design](https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0001773&viewMode=All&siblings=true) [EFO:EFO_0001773]
  - [hardware variation design](https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0001767&viewMode=All&siblings=true) [EFO:EFO_0001767]
  - [software variation design](https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0001778&viewMode=All&siblings=true) [EFO:EFO_0001778]
- [calibration](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0000818&viewMode=All&siblings=true) [OBI:OBI_0000818]
- [paired-end library preparation](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0001852&viewMode=All&siblings=true) [OBI:OBI_0001852]
  - if you did not perform paired-end library preparation, [library preparation](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0000711&viewMode=All&siblings=true) [OBI:OBI_0000711] may be more appropriate for you
- [multiplexed nucleotide library sequencing](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0001959&viewMode=All&siblings=true) [OBI:OBI_0001959]

#### geo_loc_name
Consider using the following [GAZ](https://www.ebi.ac.uk/ols/ontologies/gaz) ontology terms or terms from the [INSDC country list](http://insdc.org/country.html):
- Dalhousie: [Atlantic Ocean](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000344&lang=en&viewMode=All&siblings=false) [GAZ:GAZ_00000344];[Northwest Atlantic Ocean](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00013760&lang=en&viewMode=All&siblings=false) [GAZ:GAZ_00013760];[Bedford Basin](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00002972) [GAZ:GAZ_00002972]
- FRAM: [Atlantic Ocean](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000344&lang=en&viewMode=All&siblings=false) [GAZ:GAZ_00000344];[Greenland Sea](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00008644) [GAZ:GAZ_00008644];[Fram Strait](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000060) [GAZ:GAZ_00000060]
- MBARI: [Pacific Ocean](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000360&lang=en&viewMode=All&siblings=false) [GAZ:GAZ_00000360];[North East Pacific Ocean coastal waters of California](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00051136) [GAZ:GAZ_00051136];MBARI time series station on the shelf break
- NOAA: [Pacific Ocean](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000360&lang=en&viewMode=All&siblings=false) [GAZ:GAZ_00000360];[North East Pacific Ocean coastal waters of California](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00051136) [GAZ:GAZ_00051136];[Scripps pier](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00002545) [GAZ:GAZ_00002545]
- NOC: [Atlantic Ocean](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000344&lang=en&viewMode=All&siblings=false) [GAZ:GAZ_00000344];[English Channel coastal waters of England](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00144353) [GAZ:GAZ_00144353];[Coastal Station L4](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00187526) [GAZ:GAZ_00187526]
- SBR: [Atlantic Ocean](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000344&lang=en&viewMode=All&siblings=false) [GAZ:GAZ_00000344];[English Channel coastal waters of France](https://www.ebi.ac.uk/ols/ontologies/gaz/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00144352) [GAZ:GAZ_00144352];Roscoff-Astan observatory site/SOMLIT-Astan long-term observatory site

#### env_[broad_scale/local_scale/medium]
Consider using the following [ENVO](https://www.ebi.ac.uk/ols/ontologies/envo) terms:

Broad:
- [marine environment](https://www.ebi.ac.uk/ols/ontologies/envo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_01000320&viewMode=All&siblings=true) [ENVO:ENVO_01000320]
- [oceanic epipelagic zone biome](https://www.ebi.ac.uk/ols/ontologies/envo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_01000035&viewMode=All&siblings=false) [ENVO:ENVO_01000035]

Local:
- [sea surface layer](https://www.ebi.ac.uk/ols/ontologies/envo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_01001581) [ENVO:ENVO_01001581]
- [seasonal marine thermocline](https://www.ebi.ac.uk/ols/ontologies/envo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_01000107) [ENVO:ENVO_01000107]

Medium:
- [sea water](https://www.ebi.ac.uk/ols/ontologies/envo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_00002149) [ENVO:ENVO_00002149]
- [coastal sea water](https://www.ebi.ac.uk/ols/ontologies/envo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_00002150) [ENVO:ENVO_00002150]

#### samp_mat_processing
Consider using the following [OBI](https://www.ebi.ac.uk/ols/ontologies/obi) terms:
- [environmental material collection process](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0600012&viewMode=All&siblings=true) [OBI:OBI_0600012]
- [filtration](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0302885&viewMode=All&siblings=true) [OBI:OBI_0302885]
- [freezing storage](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0000915&viewMode=All&siblings=true) [OBI:OBI_0000915]
- [DNA extraction](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0000257&viewMode=All&siblings=true) [OBI:OBI_0000257]
- [polymerase chain reaction](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0000415&viewMode=All&siblings=true) [OBI_OBI_0000415]
- [paired-end library preparation](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0001852&viewMode=All&siblings=true) [OBI:OBI_0001852]
- [amplicon sequencing assay](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0002767&viewMode=All&siblings=true) [OBI:OBI_0002767]
- [sequence annotation](https://www.ebi.ac.uk/ols/ontologies/obi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0001944&viewMode=All&siblings=true) [OBI:OBI_0001944]

 
 
 
 
