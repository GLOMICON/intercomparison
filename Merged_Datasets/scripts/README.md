# GLOMICON scripts overview

### Order of scripts and general purpose

1. Merge individual GLOMICON datasets together: `GLOMICON_data_merge.ipynb`
	
	input files:

	> `NOAA/NOAA_GLOMICON_results/asv_taxa_sample_table.tsv`
	
	> `AWI/Glomicon-AWI-310124/GLOMICON-INTERCOMP_R-4.3.2_seqtab.merged.nochim.PR2-500.csv`
	
	> `SBR/SBR_original_files/outputs/asv_table/18s_dada2_v1.0.filtered.table.with.taxo.vsearch_BH.tsv`
	
	> `MBARI/Analysis_In_Progress/data/filtered_seq_data/GLOMICON_18S_otu_Filtered.csv`
	> `MBARI/Analysis_In_Progress/data/filtered_seq_data/GLOMICON_18S_taxa_Filtered.csv`
	> `MBARI/Analysis_In_Progress/data/filtered_seq_data/GLOMICON_18S_meta_Filtered.csv`
	> `MBARI/Analysis_In_Progress/data/filtered_seq_data/GLOMICON_18S_seq_Filtered.csv`
	
	
	> `UDalhousie/UDalhousie_GLOMICON_data/GLOMICON_UDAL_f20_feature-table_w_tax.txt`
	> `UDalhousie/UDalhousie_GLOMICON_data/GLOMICON_UDAL_f20_dna-sequences.fasta`

	output files:
	> `Merged_Datasets/data/GLOMICON_asv_merged.csv` asv table
 	> 
	> `Merged_Datasets/data/GLOMICON_taxa_merged.csv` taxonomy table
	> 
	> `Merged_Datasets/data/GLOMICON_seq_merged.csv` sequence table
	> 
	> `Merged_Datasets/data/GLOMICON_meta_merged.csv` metadata table
	> 
	> `Merged_Datasets/data/GLOMICON_merged.fasta` sequence fasta file
	

2. Run `gemelli` Robust PCA: `RPCA_Genelli.ipynb`
3. Use Dada2 and PR2 to reassign taxonomy to all ASVs:`dada2_PR2_assign_taxonomy.R`
4. Take dada2 PR2 taxonomy output, compare to blast PR2 results and choose bootstrap cutoff value, limit to eukaryotic phytoplankton, and export new taxonomy table `Format_PR2_taxonomy.ipynb` new table: `GLOMICON_asv_PR2_50_filtered.csv`
5. Parse blast xml results from blastn of ASVs against PR2, compare to genbank/MEGAN results [remove?] `make_custom_blast_taxa_table_PR2.ipynb`
6. Merge project metadata files `merge_metadata.ipynb`

7.  Plot `gemelli` program results by PC score and loading score `RPCA_bymarker.R`
8. Plot `gemelli` program results by distance metric between samples `Aitchenson_Distance_Plots.R`
9. Plot barpot overview of samples from original taxonomic labels `barplot_overview.R`
10. Plot mock community composition from original taxonomic annotations `mock_overview_orig.R` 
11. Plot mock community composition from PR2 reannotations `mock_overview_pr2.R`
12. Plot barplot overview from PR2 reannotations filtered for phytoplankton only
 `barplot_overview_pr2_filtered.R` 
13. Plot barplot overview from PR2 reannotations filtered for phytoplankton only and with replicate samples merged `barplot_replicates_merged_pr2_filtered.R`

14. `Diversity.R` Create diversity and total read plots


**Scripts/Analysis To Remove:**

1. `Check_Taxonomy.ipynb` [can remove - wanted to look at specific annotations for bugs]
2. Scripts to use taxize to convert existing taxonomy to NCBI taxonomy + plot ?

	>`Match_Taxonomy.R`
	>`Taxonomy_manual_decisions_taxizeprompts.txt`
	2. `Match_Taxonomy.R` Feed taxonomy from multiple sources through `R taxize` to put into NCBI taxonomy tree
3. Scripts to reannotate ASVs using NCBI genbank NR + MEGAN? (only use PR2?)

	>`Custom_blast_GLOMICON.sh`
	>`Format_MEGAN_output.ipynb`
	>`barplot_overview_blastnr.R`
	>`barplot_overview_blastnr_phyto.R`
	>`Limit_by_Taxonomy.ipynb`
	1. `Custom_blast_GLOMICON.sh` bash script to blastn all ASV sequences across projects to genbank NR and run through MEGAN LCA (run on virtual machine)
	2. `Format_MEGAN_output.ipynb` Takes blastn MEGAN LCA results and creates taxonomy table: `GLOMICON_taxa_blastnr.csv`
	3. `barplot_overview_blastnr.R` Plot barplot overview of community composition using blastn genbank NR results
	4. `barplot_overview_blastnr_phyto.R` Subset barplots to only phytoplankton groups from genbank NR data
	5. `Limit_by_Taxonomy.ipynb` take blastn NR taxonomy and limit to eukaryotic phytoplankton. Export: `GLOMICON_asv_limitByTaxa.csv` `GLOMICON_taxa_limitByTaxa.csv`
