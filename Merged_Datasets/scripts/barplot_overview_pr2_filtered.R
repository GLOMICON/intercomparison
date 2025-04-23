## written Katie Pitz
## 02/10/25
## Make Bar plot overview


#locations to store files:
#data files
data_directory <- './Merged_Datasets/data/pr2_filtered/'
#figures
plot_dir <- './Merged_Datasets/figures/barplot_overview_pr2_filtered/'

# Potential order to plot:
#sites <- c('EVENMOCK','BLOOMMOCK','Arctic','North Atlantic', 'English Channel', 'La Manche', 'Central California', 'Southern California')
sites <- c('EVENMOCK','BLOOMMOCK','Fram Straight','Bedford Basin', 'Western Channel', 'Roscoff', 'Monterey Bay', 'Scripps Pier')
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
library(forcats)

markers <- c("18S")
prefix = 'GLOMICON'


# Functions ---------------------------------------------------------------

make_asv_tab <- function(data_directory, marker) {
  print('ASV table')
  file = paste(prefix,"_asv_PR2_50_filtered.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  df <- read_csv(filepath) %>%
    rename_with(.cols = 1, ~"ASV")
  return(df)
}

make_taxa_tab <- function(data_directory, marker) {
  print('taxa table')
  file = paste(prefix,"_taxa_PR2_50_filtered.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  df <- read_csv(filepath) # %>%
  # rename_with(.cols = 1, ~"ASV")
  return(df)
}

make_meta_tab <- function(data_directory, marker) {
  print('metadata table')
  # file = paste(prefix,"_meta_merged.csv", sep='')
  # file = '../data/GLOMICON_meta_merged.csv'
  # filepath = paste(data_directory, file, sep='')
  filepath = './Merged_Datasets/data/GLOMICON_meta_merged.csv'
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
    facet_grid(fct_relevel(site, site_list) ~fct_relevel(Analyzing_Institute, institutes)) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    scale_x_continuous(breaks=c(1,5,10)) +
    labs(x="",y="Percent Reads of Top 20 Taxa")+
    theme_minimal() +
    theme(
      #legend
      legend.position="right",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=6,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=7,face="bold"),
      axis.text.x=element_text(size=7,colour=textcol),
      axis.text.y=element_text(size=7,colour=textcol),
      plot.background=element_blank(),
      panel.border=element_blank(),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
      plot.title=element_blank(),
      # facet_grid label text
      strip.text.y = element_text(size = 4),
      strip.text.x = element_text(size = 6))
  return(bp_top)
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
    mutate(site = case_when(Collecting_Institute == 'AWI'~ 'Fram Straight',
                            Collecting_Institute == 'MBARI'~ 'Monterey Bay',
                            Collecting_Institute == 'NOAA'~ 'Scripps Pier',
                            Collecting_Institute == 'SBR'~ 'Roscoff',
                            Collecting_Institute == 'UDalhousie'~ 'Bedford Basin',
                            Collecting_Institute == 'NOC'~ 'Western Channel',
                            TRUE ~ Collecting_Institute
    ))
  
  #OTU table long format with percent total reads
  potu.c <- make_compositional(otu.c)
  
  # Lowest Taxonomic Annotation ---------------------------------------------
  
  # species_label <- tax.c %>%
  #   unite(Class_join, Phylum, Class, sep='_', remove=FALSE) %>%
  #   unite(order_join, Phylum, Class, order, sep='_', remove=FALSE) %>%
  #   unite(family_join, Phylum, Class, order, family, sep='_', remove=FALSE) %>%
  #   unite(genus_join, Phylum, Class, order, family, genus, sep='_', remove=FALSE) %>%
  #   unite(species_join, Phylum, Class, order, family, genus, species, sep='_', remove=FALSE) %>%
  #   select(-order, -family, -genus, -species)
  
  # pr2 version:
  species_label <- tax.c %>%
    unite(class_join, division, class, sep='_', remove=FALSE) %>%
    unite(order_join, division, class, order, sep='_', remove=FALSE) %>%
    unite(family_join, division, class, order, family, sep='_', remove=FALSE) %>%
    unite(genus_join, division, class, order, family, genus, sep='_', remove=FALSE) %>%
    unite(species_join, division, class, order, family, genus, species, sep='_', remove=FALSE) %>%
    select(-order, -family, -genus, -species)
  
  
  # Overall Barplots by Taxa_level ------------------------------------------------------------
  
  ### all sites --------------
  taxas <- c('division', 'class_join', 'order_join', 'family_join', 'genus_join', 'species_join')
  for (val in taxas) {
    top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, sites)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, sites)
    barplot <- barplot_by_site(df,val, sites)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.png', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =12, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.svg', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =12, units = 'in')
    
  }
  
  ### just controls ---------------
  controls <- c('EVENMOCK','BLOOMMOCK')
  taxas <- c('division', 'class_join', 'order_join', 'family_join', 'genus_join', 'species_join')
  for (val in taxas) {
    top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, controls)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, controls)
    barplot <- barplot_by_site(df,val, controls)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_controls.png', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_controls.svg', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
  ### By Environ. group  -----------------
  #controls <- c('English Channel')
  environ_list <- c('Fram Straight','Bedford Basin', 'Western Channel', 'Roscoff', 'Monterey Bay', 'Scripps Pier')
  taxas <- c('division', 'class_join', 'order_join', 'family_join', 'genus_join', 'species_join')
  
  for (environ in environ_list) {
    controls <- c(environ)
    for (val in taxas) {
      top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, controls)
      df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, controls)
      barplot <- barplot_by_site(df,val, controls)
      #make_top20_plot(val, potu.c,species_label,meta_tab)
      
      taxa_level = sym(val)
      filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_shared_',environ,'.png', sep='')
      #print('Plot of top 20 genus average by month:')
      print(filename)
      ggsave(filename,height = 5, width =10, units = 'in')
      filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_shared',environ,'.svg', sep='')
      #print('Plot of top 20 genus average by month:')
      print(filename)
      ggsave(filename,height = 5, width =10, units = 'in')
      
    }
  }
  
  
  ### just environmental samples --------------
  controls <- c('Fram Straight','Bedford Basin', 'Western Channel', 'Roscoff', 'Monterey Bay', 'Scripps Pier')
  taxas <- c('division', 'class_join', 'order_join', 'family_join', 'genus_join', 'species_join')
  for (val in taxas) {
    top_20taxa <- make_top20_taxa(val, potu.c, species_label, meta_tab, controls)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, species_label, meta_tab, controls)
    barplot <- barplot_by_site(df,val, controls)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_environ.png', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_environ.svg', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
  
  
  # Select Specific Taxa ----------------------------------------------------
  
  ## by division -----------------
  #get top 20 most abundant taxa for dataset
  #restrict to certain Phyla for each marker?
  if (marker =='18S') {
    print(marker)
    #phyla <- c('Bacillariophyta', 'Haptophyta', 'Chlorophyta', 'unknown')
    phyla <- c('Alveolata', 'Chlorophyta', 'Cryptophyta', 'Haptophyta', 'Stramenopiles')
  }
  
  # Limit to just taxa within division
  # Plot by species
  for (taxon in phyla) {
    taxon_id = sym(taxon)
    phy_lim <- species_label %>%
      filter(division == taxon_id)
    top_20taxa <- make_top20_taxa('family_join', potu.c, phy_lim, meta_tab, sites)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, phy_lim, meta_tab, sites)
    barplot <- barplot_by_site(df,'family_join', sites)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    #taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.png', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.svg', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
  ## by class -----------------------
  
  if (marker =='18S') {
    print(marker)
    phyla <- c('Dinophyceae', 'Syndiniales','Spirotrichea', 'Prymnesiophycaea', 'Mamiellophyceae', 'Mediophyceae')
  }
  
  # Limit to just taxa within division
  # Plot by species
  for (taxon in phyla) {
    taxon_id = sym(taxon)
    phy_lim <- species_label %>%
      filter(class == taxon_id)
    top_20taxa <- make_top20_taxa('family_join', potu.c, phy_lim, meta_tab, sites)
    df <- merge_data_limit_byTopTaxa(top_20taxa, potu.c, phy_lim, meta_tab, sites)
    barplot <- barplot_by_site(df,'family_join', sites)
    #make_top20_plot(val, potu.c,species_label,meta_tab)
    
    #taxa_level = sym(val)
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.png', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxon_id,'_bar.svg', sep='')
    #print('Plot of top 20 genus average by month:')
    print(filename)
    ggsave(filename,height = 5, width =10, units = 'in')
    
  }
  
}

