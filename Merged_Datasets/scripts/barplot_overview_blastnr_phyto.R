## written Katie Pitz
## 02/10/25
## Make Bar plot overview


#locations to store files:
#data files
data_directory <- './Merged_Datasets/data/'
#figures
plot_dir <- './Merged_Datasets/figures/barplot_overview_blastnr_phyto/'

# Potential order to plot:
sites <- c('EVENMOCK','BLOOMMOCK','Arctic','North Atlantic', 'English Channel', 'La Manche', 'Central California', 'Southern California')
institutes <- c('AWI', 'SBR', 'UDAL', 'MBARI', 'NOAA')

# Load Libraries -----------------------------------------------------------------

library(readr) #read csv files
library(lubridate) #for date modifications
library(dplyr)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(tidyr)
library(stringr)
library(RColorBrewer) #colors for plotting

markers <- c("18S")
prefix = 'GLOMICON'


# Functions ---------------------------------------------------------------

make_asv_tab <- function(data_directory, marker) {
  print('ASV table')
  file = paste(prefix,"_asv_limitByTaxa.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  df <- read_csv(filepath) %>%
    rename_with(.cols = 1, ~"ASV")
  return(df)
}

make_taxa_tab <- function(data_directory, marker) {
  print('taxa table')
  file = paste(prefix,"_taxa_limitByTaxa.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  df <- read_csv(filepath) # %>%
    # rename_with(.cols = 1, ~"ASV")
  return(df)
}

make_meta_tab <- function(data_directory, marker) {
  print('metadata table')
  file = paste(prefix,"_meta_merged.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  df <- read_csv(filepath) %>% rename('SampleID' = 'sample_name') 
  return(df)
}

make_compositional <- function(df) {
  df %<>%
    tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
    group_by(SampleID) %>%
    mutate(per_tot = reads / sum(reads, na.rm=TRUE) *100) %>%
    ungroup() %>%
    arrange(-reads)
  return(df)
}

make_top20_taxa <- function(taxa_level_value, potu_df, taxa_df, meta_df, site_list) {
  taxa_level = sym(taxa_level_value)
  top_taxa_df <- potu_df %>%
    left_join(meta_df) %>%
    filter(site %in% site_list) %>%
    right_join(taxa_df) %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    filter(sum_per_tot >0) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(20)
  return(top_taxa_df)
}

merge_data_limit_byTopTaxa <- function(top_taxa_df, potu_df, taxa_df, meta_df, site_list) {
  merged_top_df <- potu_df %>% 
    full_join(taxa_df) %>% #join with taxonomy
    right_join(top_taxa_df) %>% #limit to top 20
    left_join(meta_df) %>% #join with metadata
    #filter(Collecting_Institute !='NA') %>%
    filter(site %in% site_list) %>%
    filter(replicateID <6) # don't include replicate sequenced samples (AWI)
  return(merged_top_df)
}


barplot_by_site <- function(merged_top_df,taxa_level_value, site_list) {
  taxa_level = sym(taxa_level_value)
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- merged_top_df %>%
    ggplot(aes(x = replicateID, y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = !!taxa_level))+
    facet_grid(fct_relevel(Analyzing_Institute, institutes) ~fct_relevel(site, site_list)) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    scale_x_continuous(breaks=c(1,5,10)) +
    labs(x="",y="Percent Total Reads")+
    theme_minimal() +
    theme(
      #legend
      legend.position="right",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=7,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=7,face="bold"),
      axis.text.x=element_text(size=8,colour=textcol),
      axis.text.y=element_text(size=7,colour=textcol),
      plot.background=element_blank(),
      panel.border=element_blank(),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
      plot.title=element_blank(),
      # facet_grid label text
      strip.text.x = element_text(size = 6))
  return(bp_top)
}




# plot top 20, input taxon level, potu.c,
# make_top20_plot <- function(val, potu.c,species_label,meta_tab) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    full_join(species_label) %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    filter(sum_per_tot >0) %>%
    #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    #print(n = Inf) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(20)
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- potu.c %>% 
    full_join(species_label) %>% #join with taxonomy
    right_join(top_taxa) %>% #limit to top 20
    left_join(meta_tab) %>%
    #filter(Collecting_Institute %in% c('NOC', 'UDalhousie')) %>%
    filter(Collecting_Institute !='NA') %>%
    filter(replicateID <6) %>%
    ggplot(aes(x = replicateID, y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = !!taxa_level))+
    #facet_grid(Analyzing_Institute~site) +
    facet_grid(Analyzing_Institute ~fct_relevel(site, sites)) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    scale_x_continuous(breaks=c(1,5,10)) +
    labs(x="",y="Percent Total Reads")+
    theme_minimal() +
    #facet_grid(rows = vars(PC_sign)) +
    theme(
      #legend
      legend.position="right",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=7,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=7,face="bold"),
      axis.text.x=element_text(size=8,colour=textcol),
      axis.text.y=element_text(size=7,colour=textcol),
      plot.background=element_blank(),
      panel.border=element_blank(),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
      plot.title=element_blank(),
      # facet_grid label text
      strip.text.x = element_text(size = 6))
  
  filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 5, width =10, units = 'in')
  filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.svg', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 5, width =10, units = 'in')
  
}

# plot top 20 species by Phylum
# make_top20_plot_byPhyla <- function(val, potu.c,species_label,meta_tab) {
  var = sym(val)
  taxa_level = sym('Species_join')
  top_taxa <- potu.c %>%
    full_join(species_label) %>%
    filter(Phylum == val) %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    filter(sum_per_tot >0) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(20)
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- potu.c %>% 
    full_join(species_label) %>% #join with taxonomy
    right_join(top_taxa) %>% #limit to top 20
    left_join(meta_tab) %>%
    #filter(Collecting_Institute %in% c('NOC', 'UDalhousie')) %>%
    filter(Collecting_Institute !='NA') %>%
    filter(replicateID <6) %>%
    ggplot(aes(x = replicateID, y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = !!taxa_level))+
    #facet_grid(Analyzing_Institute~site) +
    facet_grid(Analyzing_Institute ~fct_relevel(site, sites)) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    scale_x_continuous(breaks=c(1,5,10)) +
    labs(x="",y="Percent Total Reads")+
    theme_minimal() +
    #facet_grid(rows = vars(PC_sign)) +
    theme(
      #legend
      legend.position="right",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=7,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=7,face="bold"),
      axis.text.x=element_text(size=8,colour=textcol),
      axis.text.y=element_text(size=7,colour=textcol),
      plot.background=element_blank(),
      panel.border=element_blank(),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
      plot.title=element_blank(),
      # facet_grid label text
      strip.text.x = element_text(size = 6))
  
  filename = paste(plot_dir, marker,'_top20sp_bar_',var,'.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 5, width =10, units = 'in')
  filename = paste(plot_dir, marker,'_top20sp_bar_',var,'.svg', sep='')
  print(var)
  filename
  ggsave(filename,height = 5, width =10, units = 'in')
  bp_top
  
}

# Limit by Taxonomy ----------------------------------------

# Run over all markers ------------------------------------------------------

for (val in markers) {
  marker = sym(val)
  
  # Import Data -------------------------------------------------------------
  #ASV table
  print('ASV table')
  otu.c <- make_asv_tab(data_directory, marker)
  
  #taxa table
  print('taxa table')
  tax.c <- make_taxa_tab(data_directory, marker)
  
  #metadata table
  print('metadata table')
  meta_tab <- make_meta_tab(data_directory, marker)
  
  # create replicate numbers for plotting
  meta_tab %<>% 
    group_by(Analyzing_Institute, Collecting_Institute) %>%
    mutate(replicateID = row_number()) %>%
    ungroup()
  # create site names
  meta_tab %<>%
    mutate(site = case_when(Collecting_Institute == 'AWI'~ 'Arctic',
                            Collecting_Institute == 'MBARI'~ 'Central California',
                            Collecting_Institute == 'NOAA'~ 'Southern California',
                            Collecting_Institute == 'SBR'~ 'La Manche',
                            Collecting_Institute == 'UDalhousie'~ 'North Atlantic',
                            Collecting_Institute == 'NOC'~ 'English Channel',
                            TRUE ~ Collecting_Institute
    ))
  
  #OTU table long format with percent total reads
  potu.c <- make_compositional(otu.c)
  
  # Lowest Taxonomic Annotation ---------------------------------------------

  species_label <- tax.c %>%
    unite(Class_join, Phylum, Class, sep='_', remove=FALSE) %>%
    unite(Order_join, Phylum, Class, Order, sep='_', remove=FALSE) %>%
    unite(Family_join, Phylum, Class, Order, Family, sep='_', remove=FALSE) %>%
    unite(Genus_join, Phylum, Class, Order, Family, Genus, sep='_', remove=FALSE) %>%
    unite(Species_join, Phylum, Class, Order, Family, Genus, Species, sep='_', remove=FALSE) %>%
    select(-Order, -Family, -Genus, -Species)

  
  # Overall Barplots by Taxa_level ------------------------------------------------------------
  
  ### all sites --------------
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, sites)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, sites)
    barplot <- barplot_by_site(df,val, sites)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =12, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =12, units = 'in')
    
  }
  
  ### just controls ---------------
  controls <- c('EVENMOCK','BLOOMMOCK')
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, controls)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, controls)
    barplot <- barplot_by_site(df,val, controls)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_controls.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_controls.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
  ### just shared samples -----------------
  controls <- c('EVENMOCK','BLOOMMOCK', 'English Channel')
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, controls)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, controls)
    barplot <- barplot_by_site(df,val, controls)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_shared.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_shared.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
  ### just environmental samples --------------
  controls <- c('Arctic','North Atlantic', 'English Channel', 'La Manche', 'Central California', 'Southern California')
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, controls)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, controls)
    barplot <- barplot_by_site(df,val, controls)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_environ.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_environ.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
  
  
  # Select Specific Taxa ----------------------------------------------------
  
  ## by Phylum -----------------
  #get top 20 most abundant taxa for dataset
  #restrict to certain Phyla for each marker?
  if (marker =='18S') {
    print(marker)
    phyla <- c('Bacillariophyta', 'Haptophyta', 'Chlorophyta', 'unknown')
  }

  # Limit to just taxa within Phylum
  # Plot by Species
  for (taxon in phyla) {
    taxon_id = sym(taxon)
    phy_lim <- species_label %>%
      filter(Phylum == taxon_id)
    top_20taxa <- make_top20_taxa('Family_join', potu.c, phy_lim, meta_tab, sites)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, phy_lim, meta_tab, sites)
    barplot <- barplot_by_site(df,'Family_join', sites)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    #taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
  ## by Class -----------------------
  
  if (marker =='18S') {
    print(marker)
    phyla <- c('Dinophyceae', 'Cryptophyceae', 'Dictyochophyceae', 'Choanoflagellata')
  }
  
  # Limit to just taxa within Phylum
  # Plot by Species
  for (taxon in phyla) {
    taxon_id = sym(taxon)
    phy_lim <- species_label %>%
      filter(Class == taxon_id)
    top_20taxa <- make_top20_taxa('Family_join', potu.c, phy_lim, meta_tab, sites)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, phy_lim, meta_tab, sites)
    barplot <- barplot_by_site(df,'Family_join', sites)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    #taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
}

