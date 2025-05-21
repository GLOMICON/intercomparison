# GLOMICON scripts overview

### Order of scripts and general purpose

1. `GLOMICON_data_merge.ipynb` First script to merge individual GLOMICON datasets together:

	> `NOAA/NOAA_GLOMICON_results/asv_taxa_sample_table.tsv`
	
	> `AWI/Glomicon-AWI-310124/GLOMICON-INTERCOMP_R-4.3.2_seqtab.merged.nochim.PR2-500.csv`
	
	> `SBR/SBR_original_files/outputs/asv_table/18s_dada2_v1.0.filtered.table.with.taxo.vsearch_BH.tsv`
	
	> `MBARI/Analysis_In_Progress/data/filtered_seq_data/GLOMICON_18S_otu_Filtered.csv`
	> `MBARI/Analysis_In_Progress/data/filtered_seq_data/GLOMICON_18S_taxa_Filtered.csv`
	
	> `UDalhousie/UDalhousie_GLOMICON_data/GLOMICON_UDAL_f20_feature-table_w_tax.txt`

2. `Match_Taxonomy.R` Feed taxonomy from multiple sources through `taxize` to put into NCBI taxonomy tree
3. `Custom_blast_GLOMICON.sh` bash script to blastn all ASV sequences across projects to genbank NR and run through MEGAN LCA (run on virtual machine)
4. `Check_Taxonomy.ipynb` [can remove - wanted to look at specific annotations for bugs]
5. `RPCA_Genelli.ipynb` Run `gemelli` Robust PCA
5.  `RPCA_bymarker.R` Plot `gemelli` program results by PC score and loading score
6. `Aitchenson_Distance_Plots.R` Plot `gemelli` program results by distance metric between samples
7. `barplot_overview_blastnr.R` Barplot overview of community composition using blastn genbank NR results
8. `Format_MEGAN_output.ipynb`