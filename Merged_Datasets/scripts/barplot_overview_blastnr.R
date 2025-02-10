## written Katie Pitz
## 08/08/24
## Make Bar plot overview


#locations to store files:
#data files
data_directory <- './Merged_Datasets/data/'
#figures
plot_dir <- './Merged_Datasets/figures/barplot_overview/'


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
  file = paste(prefix,"_asv_merged.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  df <- read_csv(filepath) %>%
    rename_with(.cols = 1, ~"ASV")
  return(df)
}

make_taxa_tab <- function(data_directory, marker) {
  print('taxa table')
  file = paste(prefix,"_taxa_merged_updated.csv", sep='')
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
  
  
  #OTU table long format with percent total reads
  potu.c <- make_compositional(otu.c)
  
  # Lowest Taxonomic Annotation ---------------------------------------------
  
  head(tax.c)
  #get lowest taxonomic annotation level for ASVs
  # species_label <- tax.c %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='s_'| Species =='no_hit' ~as.character(Genus),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown' | Species =='g_'| Species =='no_hit'~as.character(Family),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Order),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Class),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Phylum),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character('Unknown'),
  #                              TRUE ~ as.character(Species)))
  
  species_label <- tax.c %>%
    unite(Class_join, Phylum, Class, sep='_', remove=FALSE) %>%
    unite(Order_join, Phylum, Class, Order, sep='_', remove=FALSE) %>%
    unite(Family_join, Phylum, Class, Order, Family, sep='_', remove=FALSE) %>%
    unite(Genus_join, Phylum, Class, Order, Family, Genus, sep='_', remove=FALSE) %>%
    unite(Species_join, Phylum, Class, Order, Family, Genus, Species, sep='_', remove=FALSE) %>%
    select(-Order, -Family, -Genus, -Species)

  
  # Bar by Taxa_level ------------------------------------------------------------
  #taxas <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    taxa_level = sym(val)
    top_taxa <- potu.c %>%
      full_join(species_label) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
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
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 3, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 3, width =15, units = 'in')
  }
  
  # Select Specific Taxa ----------------------------------------------------
  
  #get top 20 most abundant taxa for dataset
  #restrict to certain Phyla for each marker?
  #markers <- c("12S","COI","18S","12SnoCL", "COInoH")
  #if (marker %in% c('COI','COInoH')) {
  if (marker == 'COI') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera', 'unknown')
  }
  if (marker =='18S') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera')
  }
  if (marker =='12S') {
    print(marker)
    phyla <- c('Chordata')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    top_taxa <- potu.c %>%
      full_join(species_label) %>%
      filter(Phylum == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
      filter(sum_per_tot >0) %>%
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
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 3, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 3, width =15, units = 'in')
    bp_top
  }
  
  # Plot certain classes --------------
  
  if (marker =='18S') {
    print(marker)
    phyla <- c('Dinophyceae', 'Cryptophyceae', 'Dictyochophyceae', 'Choanoflagellata')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    top_taxa <- potu.c %>%
      full_join(species_label) %>%
      filter(Class == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
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
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 3, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 3, width =15, units = 'in')
    bp_top
  }
  
  if (marker =='18S') {
    print(marker)
    phyla <- c('Dinophyceae')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    taxa_level = sym('Genus_join')
    top_taxa <- potu.c %>%
      full_join(species_label) %>%
      filter(Class == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
      filter(sum_per_tot >0) %>%
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
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20gen_bar_',var,'.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 3, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20gen_bar_',var,'.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 3, width =15, units = 'in')
    bp_top
  }
  
}

# Just shared datasets: ----------------

for (val in markers) {
  marker = sym(val)
  
  ## Import Data -------------------------------------------------------------
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
  
  
  #OTU table long format with percent total reads
  potu.c <- make_compositional(otu.c)
  
  ## Lowest Taxonomic Annotation ---------------------------------------------
  
  head(tax.c)
  #get lowest taxonomic annotation level for ASVs
  # species_label <- tax.c %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='s_'| Species =='no_hit' ~as.character(Genus),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown' | Species =='g_'| Species =='no_hit'~as.character(Family),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Order),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Class),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Phylum),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character('Unknown'),
  #                              TRUE ~ as.character(Species)))
  
  species_label <- tax.c %>%
    unite(Class_join, Phylum, Class, sep='_', remove=FALSE) %>%
    unite(Order_join, Phylum, Class, Order, sep='_', remove=FALSE) %>%
    unite(Family_join, Phylum, Class, Order, Family, sep='_', remove=FALSE) %>%
    unite(Genus_join, Phylum, Class, Order, Family, Genus, sep='_', remove=FALSE) %>%
    unite(Species_join, Phylum, Class, Order, Family, Genus, Species, sep='_', remove=FALSE) %>%
    select(-Order, -Family, -Genus, -Species)
  
  ## Bar by Taxa_level ------------------------------------------------------------
  #taxas <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    taxa_level = sym(val)
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA')) %>%
      filter(Collecting_Institute %in% c('NOC', 'BLOOMMOCK', 'EVENMOCK')) %>%
      full_join(species_label) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
      filter(sum_per_tot >0) %>%
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
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA')) %>%
      filter(Collecting_Institute %in% c('NOC', 'BLOOMMOCK', 'EVENMOCK')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_shared.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 6, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_shared.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 6, width =10, units = 'in')
  }
  
  ## Select Specific Taxa ----------------------------------------------------
  
  #get top 20 most abundant taxa for dataset
  #restrict to certain Phyla for each marker?
  #markers <- c("12S","COI","18S","12SnoCL", "COInoH")
  #if (marker %in% c('COI','COInoH')) {
  if (marker == 'COI') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera', 'unknown')
  }
  if (marker =='18S') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA')) %>%
      full_join(species_label) %>%
      filter(Phylum == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
      filter(sum_per_tot >0) %>%
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
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA')) %>%
      filter(Collecting_Institute %in% c('NOC', 'BLOOMMOCK', 'EVENMOCK')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_shared.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_shared.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =10, units = 'in')
    bp_top
  }
  
  
}


# Only Control Samples: ---------------------------------------------------

for (val in markers) {
  marker = sym(val)
  
  ## Import Data -------------------------------------------------------------
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
  
  
  #OTU table long format with percent total reads
  potu.c <- make_compositional(otu.c)
  
  ## Lowest Taxonomic Annotation ---------------------------------------------
  
  head(tax.c)
  #get lowest taxonomic annotation level for ASVs
  # species_label <- tax.c %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='s_'| Species =='no_hit' ~as.character(Genus),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown' | Species =='g_'| Species =='no_hit'~as.character(Family),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Order),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Class),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Phylum),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character('Unknown'),
  #                              TRUE ~ as.character(Species)))
  
  species_label <- tax.c %>%
    unite(Class_join, Phylum, Class, sep='_', remove=FALSE) %>%
    unite(Order_join, Phylum, Class, Order, sep='_', remove=FALSE) %>%
    unite(Family_join, Phylum, Class, Order, Family, sep='_', remove=FALSE) %>%
    unite(Genus_join, Phylum, Class, Order, Family, Genus, sep='_', remove=FALSE) %>%
    unite(Species_join, Phylum, Class, Order, Family, Genus, Species, sep='_', remove=FALSE) %>%
    select( -Order, -Family, -Genus, -Species)
  
  ## Bar by Taxa_level ------------------------------------------------------------
  #taxas <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    taxa_level = sym(val)
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA')) %>%
      filter(Collecting_Institute %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
      full_join(species_label) %>%
      #unite(joined_taxa, !!taxa_level)
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      filter(sum_per_tot >2) %>%
      arrange(-sum_per_tot) %>%
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
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA')) %>%
      filter(Collecting_Institute %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_controls.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 6, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_controls.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 6, width =10, units = 'in')
  }
  
  ## Select Specific Taxa ----------------------------------------------------
  
  #get top 20 most abundant taxa for dataset
  #restrict to certain Phyla for each marker?
  #markers <- c("12S","COI","18S","12SnoCL", "COInoH")
  #if (marker %in% c('COI','COInoH')) {
  if (marker == 'COI') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera', 'unknown')
  }
  if (marker =='18S') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      filter(Collecting_Institute %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
      full_join(species_label) %>%
      filter(Phylum == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      filter(sum_per_tot >1) %>%
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
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA')) %>%
      filter(Collecting_Institute %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_controls.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =10, units = 'in')
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_controls.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =10, units = 'in')
    bp_top
  }
  
  
}


# Only Environmental ------------------------------------------------------

for (val in markers) {
  marker = sym(val)
  
  ## Import Data -------------------------------------------------------------
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
  
  
  #OTU table long format with percent total reads
  potu.c <- make_compositional(otu.c)
  
  ## Lowest Taxonomic Annotation ---------------------------------------------
  
  head(tax.c)
  #get lowest taxonomic annotation level for ASVs
  # species_label <- tax.c %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='s_'| Species =='no_hit' ~as.character(Genus),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown' | Species =='g_'| Species =='no_hit'~as.character(Family),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Order),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Class),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Phylum),
  #                              TRUE ~ as.character(Species))) %>%
  #   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character('Unknown'),
  #                              TRUE ~ as.character(Species)))
  
  species_label <- tax.c %>%
    unite(Class_join, Phylum, Class, sep='_', remove=FALSE) %>%
    unite(Order_join, Phylum, Class, Order, sep='_', remove=FALSE) %>%
    unite(Family_join, Phylum, Class, Order, Family, sep='_', remove=FALSE) %>%
    unite(Genus_join, Phylum, Class, Order, Family, Genus, sep='_', remove=FALSE) %>%
    unite(Species_join, Phylum, Class, Order, Family, Genus, Species, sep='_', remove=FALSE) %>%
    select(-Order, -Family, -Genus, -Species)
  
  ## Bar by Taxa_level ------------------------------------------------------------
  #taxas <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
  taxas <- c('Phylum', 'Class_join', 'Order_join', 'Family_join', 'Genus_join', 'Species_join')
  for (val in taxas) {
    taxa_level = sym(val)
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      #filter(Collecting_Institute %in% c('NOC', 'BLOOMMOCK', 'EVENMOCK')) %>%
      full_join(species_label) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
      filter(sum_per_tot >0) %>%
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
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      #filter(Collecting_Institute %in% c('NOC', 'BLOOMMOCK', 'EVENMOCK')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level),size=0, color=NA)+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_env.png', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 6, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20',taxa_level,'_bar_env.svg', sep='')
    #print('Plot of top 20 Genus average by month:')
    print(filename)
    ggsave(filename,height = 6, width =15, units = 'in')
  }
  
  ## Select Specific Taxa ----------------------------------------------------
  
  #get top 20 most abundant taxa for dataset
  #restrict to certain Phyla for each marker?
  #markers <- c("12S","COI","18S","12SnoCL", "COInoH")
  #if (marker %in% c('COI','COInoH')) {
  if (marker == 'COI') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera', 'unknown')
  }
  if (marker =='18S') {
    print(marker)
    phyla <- c('Annelida', 'Arthropoda', 'Bacillariophyta', 'Cnidaria','Rotifera')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      full_join(species_label) %>%
      filter(Phylum == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
      filter(sum_per_tot >0) %>%
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
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      #filter(Collecting_Institute %in% c('NOC', 'BLOOMMOCK', 'EVENMOCK')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_env.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_env.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =15, units = 'in')
    bp_top
  }
  # Plot certain classes --------------
  
  if (marker =='18S') {
    print(marker)
    phyla <- c('Dinophyceae', 'Cryptophyceae', 'Dictyochophyceae', 'Choanoflagellata')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      full_join(species_label) %>%
      filter(Class == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
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
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_env.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20sp_bar_',var,'_env.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =15, units = 'in')
    bp_top
  }
  
  if (marker =='18S') {
    print(marker)
    phyla <- c('Dinophyceae')
  }
  
  #plot species
  for (val in phyla) {
    var = sym(val)
    taxa_level = sym('Genus_join')
    top_taxa <- potu.c %>%
      left_join(meta_tab) %>%
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      full_join(species_label) %>%
      filter(Class == val) %>%
      # filter(!!taxa_level != 'Unknown') %>%
      # filter(!!taxa_level !='no_hit') %>%
      # filter(!!taxa_level !='unassigned') %>%
      # filter(!!taxa_level !='unknown') %>%
      # filter(!!taxa_level !='s_') %>%
      # filter(!!taxa_level !='g_') %>%
      group_by(!!taxa_level) %>%
      mutate(sum_per_tot = sum(per_tot, na.rm=TRUE)) %>%
      distinct(!!taxa_level,.keep_all = TRUE ) %>%
      arrange(-sum_per_tot) %>%
      #select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
      filter(sum_per_tot >0) %>%
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
      filter(Collecting_Institute %in% c('NOC', 'UDalhousie', 'MBARI', 'NOAA', 'AWI', 'SBR')) %>%
      #filter(Collecting_Institute %in% c('NOC', 'UDalhousie')) %>%
      ggplot(aes(x = replicateID, y = per_tot)) +
      geom_bar(stat = "identity", aes(fill = !!taxa_level))+
      facet_grid(Analyzing_Institute~Collecting_Institute) +
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
        plot.title=element_blank())
    
    filename = paste(plot_dir, marker,'_top20gen_bar_',var,'_env.png', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =15, units = 'in')
    filename = paste(plot_dir, marker,'_top20gen_bar_',var,'_env.svg', sep='')
    print(var)
    filename
    ggsave(filename,height = 6, width =15, units = 'in')
    bp_top
  }
  
}

