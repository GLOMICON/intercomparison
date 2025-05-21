## written Katie Pitz
## 02/10/25
## Make Bar plot overview


#locations to store files:
#data files
data_directory <- './Merged_Datasets/data/pr2_filtered/'
#figures
plot_dir <- './Merged_Datasets/figures/mock_overview_pr2/'

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
marker <- '18S'

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


# Import Data ----------------------------------------

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

# pr2 version:
species_label <- tax.c %>%
  unite(class_join, division, class, sep='_', remove=FALSE) %>%
  unite(order_join, division, class, order, sep='_', remove=FALSE) %>%
  unite(family_join, division, class, order, family, sep='_', remove=FALSE) %>%
  unite(genus_join, division, class, order, family, genus, sep='_', remove=FALSE) %>%
  unite(species_join, division, class, order, family, genus, species, sep='_', remove=FALSE) %>%
  select(-order, -family, -genus, -species)

# Handle Mock Community -----------------------

# import mock community composition
filepath = "/Users/kpitz/github/GLOMICON/intercomparison/Merged_Datasets/data/mock_composition/AWI_mock_community_composition_toPR2.csv"
print(filepath)
mockcc <- read_csv(filepath) %>%
  rename(mock_species = `Total volume ~42 µl`) %>%
  select(-`Total volume ~92 µl`)

# taxa table
tax_mockcc <- mockcc %>%
  select(mock_species,domain, supergroup, division, subdivision, class, order, family, genus, species)
# potu table plus basic metadata to match other samples
otu_mockcc <- mockcc %>%
  select(mock_species, EMper_total_copy_numbers, BMper_total_copy_numbers) %>%
  rename(EM = EMper_total_copy_numbers, BM =BMper_total_copy_numbers ) %>%
  pivot_longer(!mock_species,names_to = "SampleID", values_to = "per_tot") %>%
  mutate(Analyzing_Institute = 'Original_Concentration') %>%
  mutate(site = case_when(SampleID == 'BM' ~ 'BLOOMMOCK',
                          SampleID == 'EM' ~ 'EVENMOCK')) %>%
  mutate(replicateID = 1)

# Link mock cc to sequenced mocks -------------------

# potu, taxa table limited to mocks
df<- potu.c %>%
  left_join(meta_tab) %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
  left_join(tax.c) # %>%
  # filter(per_tot >0)
  
# plot at taxonomic levels:

site_list <- c('EVENMOCK', 'BLOOMMOCK')
institutes <- c('AWI', 'SBR', 'UDAL', 'MBARI', 'NOAA', 'Original_Concentration')
textcol <- 'grey'

taxas <- c('division', 'subdivision','class', 'order')
for (val in taxas) {
  tax_level = sym(val)
  bplot <- otu_mockcc %>%
    full_join(tax_mockcc) %>%
    bind_rows(df) %>%
    # group by tax_level in order to limit by overall abundance:
    group_by(!!tax_level, SampleID, replicateID, site, Analyzing_Institute) %>%
    mutate(per_tot = sum(per_tot)) %>%
    ungroup() %>%
    distinct(!!tax_level, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
    #total sum per_tot
    group_by(!!tax_level) %>%
    mutate(overall_abund = sum(per_tot)) %>%
    ungroup() %>%
    filter(overall_abund >3) %>%
    ggplot(aes(x = replicateID, y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = !!tax_level))+
    facet_grid(fct_relevel(site, site_list) ~fct_relevel(Analyzing_Institute, institutes)) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    scale_x_continuous(breaks=c(1,5,10)) +
    labs(x="",y="Percent of Total Reads per Sample")+
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
  
  bplot
  
  filename = paste(plot_dir, marker,'_mockcc_',tax_level,'_bar_controls.png', sep='')
  #print('Plot of top 20 genus average by month:')
  print(filename)
  ggsave(filename,height = 5, width =10, units = 'in')
}

# look at species within divisions
taxas <- c('Alveolata', 'Haptophyta','Chlorophyta', 'Prasinodermophyta', 'Stramenopiles')
for (val in taxas) {
  print(val)
  #val = 'Haptophyta'
  tax_level = sym(val)
  bplot <- otu_mockcc %>%
    full_join(tax_mockcc) %>%
    bind_rows(df) %>%
    filter(division==val) %>%
    #glimpse() %>%
    # join several tax rows:
    unite(label, division, class, family,genus, species, sep = '_') %>%
    # group by tax_level in order to limit by overall abundance:
    group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
    mutate(per_tot = sum(per_tot)) %>%
    ungroup() %>%
    distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
    #total sum per_tot
    group_by(label) %>%
    mutate(overall_abund = sum(per_tot)) %>%
    ungroup() %>%
    filter(overall_abund >1) %>%
    ggplot(aes(x = replicateID, y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = label))+
    facet_grid(fct_relevel(site, site_list) ~fct_relevel(Analyzing_Institute, institutes)) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    scale_x_continuous(breaks=c(1,5,10)) +
    labs(x="",y="Percent of Total Reads per Sample", fill=val)+
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
  
  bplot
  
  filename = paste(plot_dir, marker,'_',val,'_mockcc_bar_sp.png', sep='')
  #print('Plot of top 20 genus average by month:')
  print(filename)
  ggsave(filename,height = 5, width =10, units = 'in')

}

val = 'Stramenopiles'
print(val)
#val = 'Haptophyta'
tax_level = sym(val)
bplot <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  filter(division==val) %>%
  #glimpse() %>%
  # join several tax rows:
  unite(label, division, class, family,genus ,sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >1) %>%
  ggplot(aes(x = replicateID, y = per_tot)) +
  geom_bar(stat = "identity", aes(fill = label))+
  facet_grid(fct_relevel(site, site_list) ~fct_relevel(Analyzing_Institute, institutes)) +
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  #scale_fill_brewer(palette = "Greens", type = c("regular"), direction = 1)+
  scale_x_continuous(breaks=c(1,5,10)) +
  labs(x="",y="Percent of Total Reads per Sample", fill=val)+
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

bplot

filename = paste(plot_dir, marker,'_',val,'_mockcc_bar_gen.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 5, width =10, units = 'in')

# outside of forloop:
# join taxonomic columns to get context
# look at genus level?
#val = 'Haptophyta'
#tax_level = sym(val)
bplot <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class, order, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  ggplot(aes(x = replicateID, y = per_tot)) +
  geom_bar(stat = "identity", aes(fill = label))+
  facet_grid(fct_relevel(site, site_list) ~fct_relevel(Analyzing_Institute, institutes)) +
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  scale_x_continuous(breaks=c(1,5,10)) +
  labs(x="",y="Percent of Total Reads per Sample", fill='division_class_order')+
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

bplot

filename = paste(plot_dir, marker,'_order_mockcc_bar_3per.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 5, width =10, units = 'in')

# line plot ------------
bplot <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  filter(site=='BLOOMMOCK') %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  #ggplot(aes(x = Analyzing_Institute, y = per_tot, color=label, fill=label)) +
  ggplot(aes(x = label, y = per_tot, color=Analyzing_Institute)) +
  #geom_point() +
  geom_boxplot()

bplot


test <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  filter(site=='BLOOMMOCK') %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class,order, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  # need to pull out original conc:
  mutate(orig_conc = case_when(Analyzing_Institute=='Original_Concentration' ~ per_tot)) %>%
  filter(orig_conc >0) %>%
  select(label, orig_conc)
test

# scatter
bplot <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  filter(site=='BLOOMMOCK') %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  full_join(test) %>%
  # need to pull out original conc:
  #mutate(orig_conc = case_when(Analyzing_Institute=='Original_Concentration' ~ per_tot)) %>%
  #ggplot(aes(x = Analyzing_Institute, y = per_tot, color=label, fill=label)) +
  ggplot(aes(y = per_tot, x = orig_conc, color=Analyzing_Institute)) +
  geom_line(data=. %>% filter(Analyzing_Institute=='Original_Concentration')) +
  geom_point(aes(shape=Analyzing_Institute), alpha=0.8)

bplot

# difference from original
test <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  #filter(site=='BLOOMMOCK') %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class,order, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  # need to pull out original concs:
  filter(Analyzing_Institute == 'Original_Concentration') %>%
  mutate(BMorig_conc = case_when(SampleID=='BM' ~ per_tot)) %>%
  mutate(EMorig_conc = case_when(SampleID=='EM' ~ per_tot)) %>%
  mutate(orig_conc = case_when(SampleID %in% c('BM', 'EM') ~ per_tot)) 
  #filter(orig_conc >0) %>%
  # select(label, orig_conc)
test


means <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  #filter(site=='BLOOMMOCK') %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class,order, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  full_join(test %>% select(label, site,orig_conc)) %>%
  # get difference from original:
  mutate(diff = per_tot - orig_conc) %>%
  filter(orig_conc >0) %>%
  group_by(label, Analyzing_Institute, site) %>%
  mutate(mean = mean(diff)) %>%
  ungroup() %>%
  distinct(label, Analyzing_Institute, site,mean)



bplot <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  #filter(site=='BLOOMMOCK') %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class,order, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  full_join(test %>% select(label, site,orig_conc)) %>%
  # get difference from original:
  mutate(diff = per_tot - orig_conc) %>%
  filter(orig_conc >0) %>%
  # need to pull out original conc:
  #mutate(orig_conc = case_when(Analyzing_Institute=='Original_Concentration' ~ per_tot)) %>%
  #ggplot(aes(x = Analyzing_Institute, y = per_tot, color=label, fill=label)) +
  #ggplot(aes(y = diff, x = label, color=Analyzing_Institute)) +
  ggplot(aes(x = diff, y = label, color=Analyzing_Institute)) +
  #geom_line(data=. %>% filter(Analyzing_Institute=='Original_Concentration')) +
  #geom_jitter(aes(shape=Analyzing_Institute), alpha=0.8) +
  #geom_boxplot() + 
  geom_jitter(size=1) +
  geom_point(data=means, aes(x=mean, y=label, color=Analyzing_Institute, shape=Analyzing_Institute), alpha=0.6, size=4) +
  facet_grid(~site) +
  theme(axis.text.y = element_text(size=5))

bplot
filename = paste(plot_dir, marker,'_order_mockcc_scatter_horiz.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6, width =12, units = 'in')

# bar version
#palette_Bar <- c(wes_palette("Chevalier1", type = "discrete")[2], wes_palette("Darjeeling2", type = "discrete")[2], 'darkgrey')
library(wesanderson)
palette_Bar <- c(wes_palette("Darjeeling1", type = "discrete"), 'darkgrey')

bplot <- otu_mockcc %>%
  full_join(tax_mockcc) %>%
  bind_rows(df) %>%
  #filter(site=='BLOOMMOCK') %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
  #filter(division==val) %>%
  # join several tax rows:
  unite(label, division, class,order, sep = '_') %>%
  # group by tax_level in order to limit by overall abundance:
  group_by(label, SampleID, replicateID, site, Analyzing_Institute) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(label, SampleID, replicateID, site, Analyzing_Institute, per_tot) %>%
  #total sum per_tot
  group_by(label) %>%
  mutate(overall_abund = sum(per_tot)) %>%
  ungroup() %>%
  filter(overall_abund >3) %>%
  full_join(test %>% select(label, site,orig_conc)) %>%
  # get difference from original:
  mutate(diff = per_tot - orig_conc) %>%
  filter(orig_conc >0) %>%
  # need to pull out original conc:
  #mutate(orig_conc = case_when(Analyzing_Institute=='Original_Concentration' ~ per_tot)) %>%
  #ggplot(aes(x = Analyzing_Institute, y = per_tot, color=label, fill=label)) +
  #ggplot(aes(y = diff, x = label, color=Analyzing_Institute)) +
  filter(Analyzing_Institute!='Original_Concentration') %>%
  ggplot(aes(x = diff, y = label, color=fct_relevel(Analyzing_Institute, institutes))) +
  #geom_line(data=. %>% filter(Analyzing_Institute=='Original_Concentration')) +
  #geom_jitter(aes(shape=Analyzing_Institute), alpha=0.8) +
  #geom_boxplot() + 
  geom_jitter(size=0.5) +
  #geom_point(data=means, aes(x=mean, y=label, color=Analyzing_Institute, shape=Analyzing_Institute), alpha=0.6, size=4) +
  geom_bar(data=means %>% filter(Analyzing_Institute!='Original_Concentration'), stat='identity',position=position_dodge(width=0.95), aes(x=mean, y=label, fill=fct_relevel(Analyzing_Institute, institutes)), alpha=0.8, width=0.8) +
  facet_grid(~site) +
  #scale_fill_brewer(palette = "Dark2", type = c("regular"), direction = 1)+
  #scale_color_brewer(palette = "Dark2", type = c("regular"), direction = 1)+
  scale_fill_manual(values = palette_Bar) +
  scale_color_manual(values = palette_Bar) +
  #facet_grid(site~Analyzing_Institute) +
  labs(x='Difference (Measured - True) Mock Percent Abundance', fill='Analyzing Institute', color='Analyzing Institute') +
  theme(axis.text.y = element_text(size=5))

bplot
filename = paste(plot_dir, marker,'_order_mockcc_bar_horiz.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =12, units = 'in')


# Radar plots ----------
library(fsmb)

# Parallel coordinates ---------
library(GGally)
test <- mockcc %>%
  rename(EM = EMper_total_copy_numbers, BM = BMper_total_copy_numbers) %>%
  select(mock_species, EM, BM) %>%
  ggparcoord(groupColumn=1, columns=2:3,showPoints = TRUE )
test


#ggparcoord()

# SCRAP below ---------------


# # get list of genera:
# mock_genera <- mockcc %>%
#   group_by(genus) %>%
#   mutate(BM = sum(BMper_total_copy_numbers)) %>%
#   mutate(EM = sum(EMper_total_copy_numbers)) %>%
#   ungroup() %>%
#   distinct(genus, BM, EM, .keep_all = TRUE) %>%
#   # remove species column though:
#   select(!species) %>%
#   pivot_longer(!genus,names_to = "SampleID", values_to = "per_tot") %>%
#   mutate(Analyzing_Institute = 'Original_Concentration') %>%
#   mutate(site = case_when(SampleID == 'BM' ~ 'BLOOMMOCK',
#                           SampleID == 'EM' ~ 'EVENMOCK')) %>%
#   mutate(replicateID = 1)
# mock_genera

# Select only taxa which match mock communities
site_list <- c('EVENMOCK', 'BLOOMMOCK')
institutes <- c('AWI', 'SBR', 'UDAL', 'MBARI', 'NOAA', 'Original_Concentration')
textcol <- 'grey'
marker <- '18S'

#otu taxa table for mock cc:
mock_potu <- tax.c %>%
  right_join(mock_genera) %>%
  distinct(genus, SampleID, .keep_all = TRUE)
# otu table limited by top genera
df<- tax.c %>%
  right_join(mock_genera %>% select(genus)) %>%
  # select(-SampleID, -per_tot) %>%
  left_join(potu.c) %>%
  #full_join(mock_potu) %>%
  left_join(meta_tab) %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK'))
# now join both together and plot:
bplot <- df %>%
  full_join(mock_potu) %>%
  ggplot(aes(x = replicateID, y = per_tot)) +
  geom_bar(stat = "identity", aes(fill = genus))+
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
bplot






test <- tax.c %>%
  right_join(mock_genera) %>%
  left_join(potu.c) %>%
  left_join(meta_tab) %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK', 'BM', 'EM')) %>%
  ggplot(aes(x = replicateID, y = per_tot)) +
  geom_bar(stat = "identity", aes(fill = genus))+
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
  
test






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

