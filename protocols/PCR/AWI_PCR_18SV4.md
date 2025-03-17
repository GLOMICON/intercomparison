---
# MIOP terms
methodology_category: Omics Analysis
project: # can be general
purpose: PCR [OBI:0000415]
analyses: PCR [OBI:0000415]
geographic_location: Arctic Ocean [GAZ:00000323]
broad_scale_environmental_context: marine biome [ENVO:00000447], marine photic zone [ENVO:00000209]
local_environmental_context: marine biome [ENVO:00000447], marine photic zone [ENVO:00000209]
environmental_medium: sea water [ENVO:00002149], polymerase chain reaction [OBI:0000415]
target: 18S [NCIT:C48172]
creator: Christian Wolf, Felix Janssen
materials_required: vortexer [OBI:0400118], PCR instrument [OBI:0000989], agarose gel electrophoresis system [OBI:0001134]
skills_required: sterile technique, pipetting skills, standard molecular technique
time_required: # minutes (integer)
personnel_required: 1
language: en
issued: # YYYY-MM-DD
audience: scientists
publisher: Alfred Wegener Institute
hasVersion: # 1
license: # CC0 1.0 Universal
maturity level: mature

# FAIRe terms
pcr_0_1: 1 
thermocycler: Eppendorf, Mastercycler # name of thermocycler
amplificationReactionVolume: 25
assay_name: # ssu16sv4v5_emp
assay_validation: # not provided
targetTaxonomicAssay: 18S rRNA gene sequencing targeting the V4 region using primers 528F and 964iR
targetTaxonomicScope: Protists
target_gene: 18S rRNA
target_subfragment: V4
ampliconSize: # 411
pcr_primer_forward: GCGGTAATTCCAGCTCCAA
pcr_primer_reverse: ACTTTCGTTCTTGATYRR
pcr_primer_name_forward: 528F
pcr_primer_name_reverse: 964iR
pcr_primer_reference_forward: 10.1093/oxfordjournals.molbev.a040362
pcr_primer_reference_reverse: not applicable
pcr_primer_vol_forward: 5
pcr_primer_vol_reverse: 5
pcr_primer_conc_forward: 1
pcr_primer_conc_reverse: 1
probeReporter: not applicable
probeQuencher: not applicable
probe_seq: not applicable
probe_ref: not applicable
probe_conc: not applicable
commercial_mm: 2x KAPA HiFi HotStart ReadyMix (Roche) # AmpliTaq Gold 360 Master Mix
custom_mm: # PCR reactions were run in 25 uL reaction volumes, with 1.0 uL of DNA, 12.5 uL of AmpliTaq Gold, 9.5 uL of water, and 1.0 uL of each primer (10 uM)
pcr_dna_vol: 2.5
pcr_rep: # 1
nucl_acid_amp: # https://doi.org/10.1111/1462-2920.13023
pcr_cond: initial denaturation:95_3;denaturation:95_0.5;annealing:55_0.5;elongation:72_0.5;final elongation:72_5;25
annealingTemp: 55
pcr_cycles: 25
pcr_analysis_software: # not provided
pcr_method_additional: # not provided
---

# AWI PCR Protocol 18S rRNA V4

## PROTOCOL INFORMATION

### Minimum Information about an Omics Protocol (MIOP)

- MIOP terms are listed in the YAML frontmatter of this page.
- See [MIOP_definition.md](https://github.com/BeBOP-OBON/0_protocol_collection_template/blob/main/MIOP_definition.md) for list and definitions.

### Making eDNA FAIR (FAIRe)

- FAIRe terms are listed in the YAML frontmatter of this page.
- See <https://fair-edna.github.io/download.html> for the FAIRe checklist and more information.
- See <https://fair-edna.github.io/guidelines.html#missing-values> for guidelines on missing values that can be used for missing FAIRe or MIOP terms.

### Authors

- All authors known to have contributed to the preparation of this protocol, including those who filled in the template.
- Visit https://orcid.org/ to register for an ORCID.
- Date is the date the author first worked on the protocol.

| PREPARED BY | AFFILIATION | ORCID | DATE |
| ------------- | ------------- | ------------- | ------------- |
| Christian Wolf  | AWI  | Content Cell | yyyy-mm-dd |
| Felix Janssen  | AWI  | Content Cell | yyyy-mm-dd |

### Related Protocols

- This section contains protocols that should be known to users of this protocol.
- Include the link to each protocol.
- Include the version number and release date (if available).
- Internal/External: "Internal" are derivative or altered protocols, or other protocols in this workflow. "External" are protcols from manufacturers or other groups.

| PROTOCOL NAME | LINK         | VERSION      | RELEASE DATE | INTERNAL/EXTERNAL |
| ------------- | ------------ | ------------ | ------------ | ----------------- |
| Place_extraction_protocol_here  | Content Cell | Content Cell | yyyy-mm-dd   | Content Cell      |
| Place_sequencing_protocol_here  | Content Cell | Content Cell | yyyy-mm-dd   | Content Cell      |

### Protocol Revision Record

- Version numbers start at 1.0.0 when the protocol is first completed and will increase when changes that impact the outcome of the procedure are made (patches: 1.0.1; minor changes: 1.1.0; major changes: 2.0.0).
- Release date is the date when a given protocol version was finalised.
- Description of revisions includes a brief description of what was changed relative to the previous version.

| VERSION | RELEASE DATE | DESCRIPTION OF REVISIONS |
| ------------- | ------------- | ------------- |
| 1.0.0 | yyyy-mm-dd | Initial release |

### Acronyms and Abbreviations

| ACRONYM / ABBREVIATION | DEFINITION |
| ------------- | ------------- |
| AWI  | Alfred Wegener Institute  |
| Content Cell  | Content Cell  |

### Glossary

| SPECIALISED TERM | DEFINITION |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |

## BACKGROUND

This document describes the required protocol to conduct [insert name of the method/protocol].

### Summary

Insert a short description of the background for the method/protocol (e.g. why and for which purpose do you perform water sampling).
Please provide a brief summary of your method including, as appropriate, a brief description of what techniques your best practice is about, which ocean environments or regions it targets, the primary sensors covered, what type of data/measurements/observing platform it covers, limits to its applicability.

### Method Description and Rationale

Insert a short description of the functioning principal of the methodology used in the protocol (i.e. how does the method work?). Please note that this is different from the step-by-step description of the protocol procedure.
Insert a short statement explaining why the specific methodology used in the protocol has been selected (e.g. it is highly reproducible, highly accurate, procedures are easy to execute etc….).

### Spatial Coverage and Environment(s) of Relevance

If applicable, please specify the region where the protocol is applied. For regional term guidance see the [GAZ ontology](https://www.ebi.ac.uk/ols4/ontologies/gaz). If applicable, please indicate here the environment(s) of relevance for the protocol, e.g. Abyssal plain. Select from the [ENVO ontology](https://www.ebi.ac.uk/ols4/ontologies/envo).

## PERSONNEL REQUIRED

Insert the number of technicians, data managers, and scientists required for the good execution of the procedure

### Safety

Identify hazards associated with the procedure and specify protective equipment and safety training required to safely execute the procedure

### Training Requirements

Specify technical training required for the good execution of the procedure.

### Time Needed to Execute the Procedure

Specify how much time is necessary to execute the procedure.

## EQUIPMENT

- Opentrons Consumables: If using Opentrons OT-2 Robot for KF Plate Prep.
- Description: E.g., "filter".
- Product Name and Model: Provide the official name of the product.
- Manufacturer: Provide the name of the manufacturer of the product.
- Quantity: Provide quantities necessary for one application of the standard operating procedure (e.g., number of filters).
- Remark: For example, some of the consumable may need to be sterilized, some commercial solution may need to be diluted or shielded from light during the operating procedure.

| DESCRIPTION | PRODUCT NAME AND MODEL | MANUFACTURER | QUANTITY | REMARK |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| **Durable equipment** |
| Thermal cycler | Mastercycler | Eppendorf | 1 | Can be substituted with generic |
| Content Cell | Content Cell | Content Cell | Content Cell | Content Cell |
| **Consumable equipment** |
| PCR Master Mix | 2x KAPA HiFi HotStart ReadyMix | Roche | 12.5 | ul |
| Content Cell | Content Cell | Content Cell | Content Cell | Content Cell |
| Content Cell | Content Cell | Content Cell | Content Cell | Content Cell |
| **Chemicals** |
| Content Cell | Content Cell | Content Cell | Content Cell | Content Cell |
| Content Cell | Content Cell | Content Cell | Content Cell | Content Cell |

## STANDARD OPERATING PROCEDURE

In the following SOP, please use the exact names of equipment as noted in the table above.

Provide a step-by-step description of the protocol. The identification of difficult steps in the protocol and the provision of recommendations for the execution of those steps are encouraged.

### Preparation

Please specify the preparatory actions you took before you collected the samples and note what equipment was needed to do so (e.g. disinfection of work surfaces, preparations to the equipment you intend to use later on).

1. [Step 1]
2. [Step 2]

### PCR

1. Make PCR master mix and add 22.5 ul to each well of PCR plate.


**Primers**: PCR primer sequences

| PCR Primer Name | Direction | Sequence (5’ -> 3’)|
| ----- | ----- | ----- |
| 528F | forward | GCGGTAATTCCAGCTCCAA |
| 964iR | reverse | ACTTTCGTTCTTGATYRR |

**Reaction Mixture**: PCR reagents, volumes, initial and final concentrations

| Reagent | Volume | Initial Concentration | final concentration|
| ----- | ----- | ----- | ----- |
| 2x KAPA HiFi HotStart ReadyMix (Roche) | 12.5ul | not applicable |not applicable |
| Fwd primer | 5 ul | 1 μM | 0.2 μM |
| Rev primer | 5 ul | 1 μM | 0.2 μM |

2. Add 2.5 ul of template DNA to respective wells for a total reaction volume of 25 ul per well.
3. [next step]
4. [next step]


**PCR Cycling Program**: 

| PCR Step | Temperature | Duration | Repetition |
| ----- | ----- | ----- | ----- |
| Initial Denaturation | 95°C | 3min | 1x |
| Denaturation | 95°C | 30s | 25x |
| Annealing | 55°C | 30s | 25x |
| Extension | 72°C | 30s | 25x |
| Final Extension | 72°C | 5min | 1x |
| Hold | 4°C | ∞ | |


### Quality Control

Please specify the actions you took to confirm the quality of the PCR output, to clean up the PCR output and the equipment you used (e.g. agarose gel to confirm quality, purification of PCR products).

1. [Step 1]
2. [Step 2]

#### Positive Control

Please include information about any positive controls, used in every PCR run to verify success of the PCR reaction. This should include a description of the sequence(s), the concentration and volume added, and the reference sequence(s).

#### Negative Control

Please include information about any negative controls, such as PCR-grade water used as a no template control (NTC) when setting up each PCR plate.

### PCR Clean-up

Please specify the actions you took to clean up the PCR output and which equipment you used for this (e.g. agarose gel to confirm quality, purification of PCR products).

1. The amplicons were purified away from free primers and primer dimer species with AMPure XP beads (Beckmann Coulter).
2. [Step 2]

### Basic Troubleshooting Guide

- Identify known issues associated with the procedure, if any.
- Provide troubleshooting guidelines when available.

## REFERENCES

- Insert all references cited in the document.
- Please insert full DOI address when available, e.g. http://doi.dx.org/10.1007/s11258-014-0404-1.

## APPENDIX A: DATASHEETS

Link templates (e.g. preformatted spreadsheets) used to record measurements and report on the quality of the data as well as any documents such as manufacturer specifications, images, etc that support this protocol. Please include a short note describing the document's relevance.
