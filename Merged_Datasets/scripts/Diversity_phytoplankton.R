# Proportion of reads annotated
## kpitz
## 04/22/25


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

# Set Constants -------------------------------------------------

marker <- '18S'
prefix <- 'GLOMICON'

#locations to store files:
#data files
data_directory <- './Merged_Datasets/data/'
#figures
plot_dir <- './Merged_Datasets/figures/diversity_phytoplankton/'

# Order to plot sites and institutes:
site_list <- c('EVENMOCK','BLOOMMOCK','Fram Straight','Bedford Basin', 
           'Western Channel', 'Roscoff', 'Monterey Bay', 'Scripps Pier')
#institutes <- c('AWI', 'SBR', 'UDAL', 'MBARI', 'NOAA')
institutes <- c('AWI', 'SBR', 'UDAL', 'MBARI', 'AOML', 'PARADA')

# Functions ---------------------------------------------------------------

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
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

# Import Data -------------------------------------------------------------

# asv table
#filepath = paste(data_directory, prefix, '_asv_merged.csv', sep='')
filepath = paste(data_directory,'pr2_reassigned/phytoplankton_only/',prefix, '_asv_PR2_50_filtered.csv', sep='')
otu.c <- read_csv(filepath) %>%
  rename_with(.cols = 1, ~"ASV")

# metadata table
#filepath = paste(data_directory,'merged_original_taxonomy/', prefix, '_meta_merged.csv', sep='')
filepath = paste(data_directory, prefix, '_meta_merged.csv', sep='')

meta_tab <- read_csv(filepath) %>%
  rename('SampleID' = 'sample_name') %>%
  # remove duplicate sequenced AWI samples with lower reads/diversity
  # Framstrait_01_0008 Framstrait_01_0049
  mutate(plate_string = 'none') %>%
  mutate(plate_string = case_when(Analyzing_Institute == 'AWI' ~ str_match(SampleID, '.*_.*_(.*)')[,2], TRUE ~ plate_string)) %>%
  filter(plate_string !='0008') %>%
  group_by(Analyzing_Institute, Collecting_Institute) %>%
  mutate(replicateID = row_number()) %>%
  ungroup() %>%
  mutate(site = case_when(Collecting_Institute == 'AWI'~ 'Fram Straight',
                          Collecting_Institute == 'MBARI'~ 'Monterey Bay',
                          Collecting_Institute == 'NOAA'~ 'Scripps Pier',
                          Collecting_Institute == 'AOML'~ 'Scripps Pier',
                          Collecting_Institute == 'SBR'~ 'Roscoff',
                          Collecting_Institute == 'UDalhousie'~ 'Bedford Basin',
                          Collecting_Institute == 'UDAL'~ 'Bedford Basin',
                          Collecting_Institute == 'NOC'~ 'Western Channel',
                          TRUE ~ Collecting_Institute
  )) %>%
  # remove AWI's blank sample
  filter(site !='NA') %>%
  mutate(Analyzing_Institute = case_when(Analyzing_Institute == 'NOAA' ~ 'AOML',
                                         TRUE ~ Analyzing_Institute))


# pr2 taxonomy - unfiltered (includes metazoa and fungi)
filepath = paste(data_directory, 'pr2_reassigned/phytoplankton_only/',prefix, '_taxa_PR2_50_filtered.csv', sep='')
tax.c <-read_csv(filepath)

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

# Plot total reads per sample ------------------------------------------

p <- potu.c %>%
  filter(reads>0) %>%
  group_by(SampleID) %>%
  mutate(reads = sum(reads, na.rm=TRUE)) %>%
  ungroup() %>%
  distinct(SampleID, reads) %>%
  right_join(meta_tab) %>%
  #filter(site !='NA') %>%
  #filter(plate_string !='0008') %>%
  #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
  ggplot(aes(y=reads, x=fct_relevel(site, site_list), color=fct_relevel(site, site_list), group=fct_relevel(site, site_list))) +
  geom_boxplot() +
  geom_jitter() +
  scale_color_brewer(palette = 'Dark2') +
  # facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  facet_grid(fct_relevel(Analyzing_Institute, institutes)~., scales="free") +
  labs(x='site', y='reads', color='site') +
  ggtitle('Reads per Sample') +
  theme(axis.text.x = element_text(angle=90, vjust= 0.5))

p

filename = paste(plot_dir, marker,'_totreads_byproject.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =8, units = 'in')


p <- potu.c %>%
  filter(reads>0) %>%
  group_by(SampleID) %>%
  mutate(reads = sum(reads, na.rm=TRUE)) %>%
  ungroup() %>%
  distinct(SampleID, reads) %>%
  right_join(meta_tab) %>%
  #filter(site !='NA') %>%
  #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
  ggplot(aes(y=reads, x=fct_relevel(site, site_list), color=fct_relevel(Analyzing_Institute, institutes))) +
  geom_boxplot() +
  geom_jitter() +
  scale_color_brewer(palette = 'Set1') +
  # facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  #facet_grid(fct_relevel(Analyzing_Institute, institutes)~.) +
  labs(x='site', y='reads', color='Institute') +
  ggtitle('Reads per Sample') +
  theme(axis.text.x = element_text(angle=90, vjust= 0.5))

p

filename = paste(plot_dir, marker,'_totreads_byproject_sharedy.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

# # free y axis
# p <- potu.c %>%
#   group_by(SampleID) %>%
#   mutate(reads = sum(reads, na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(SampleID, reads) %>%
#   left_join(meta_tab) %>%
#   #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
#   ggplot(aes(y=reads, x='', color=site)) +
#   geom_point() +
#   #geom_boxplot() +
#   facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
# p
# 
# 
# 
# # same y-axis
# p <- potu.c %>%
#   group_by(SampleID) %>%
#   mutate(reads = sum(reads, na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(SampleID, reads) %>%
#   left_join(meta_tab) %>%
#   #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
#   ggplot(aes(y=reads, x='', color=site)) +
#   geom_point() +
#   #geom_boxplot() +
#   facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list))
# p
# 
# # look at samples with <=5k reads
# p <- potu.c %>%
#   group_by(SampleID) %>%
#   mutate(reads = sum(reads, na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(SampleID, reads) %>%
#   # look at where we have samples with <=5k reads
#   filter(reads<=5000) %>%
#   left_join(meta_tab) %>%
#   #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
#   ggplot(aes(y=reads, x='', color=site)) +
#   geom_point() +
#   #geom_boxplot() +
#   facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list))
# p

# Plot total number of ASVs ------------------------------------------

# free y axis
p <- potu.c %>%
  # remove any nan or 0 read values
  filter(reads>0) %>%
  mutate(count = 1) %>%
  group_by(SampleID) %>%
  mutate(count = sum(count)) %>%
  ungroup() %>%
  distinct(SampleID, count) %>%
  right_join(meta_tab) %>%
  filter(site !='NA') %>%
  #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
  ggplot(aes(y=count, x=fct_relevel(site, site_list), color=fct_relevel(site, site_list))) +
  geom_boxplot() +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') +
  #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
  facet_wrap(fct_relevel(Analyzing_Institute, institutes)~.,  ncol=2, scales="free") +
  #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  labs(y='Number of ASVs', x='', color='site')+
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))
  #theme(axis.text.x = element_text(angle=45, vjust= 0.3, hjust = 0.5))
p

filename = paste(plot_dir, marker,'_totASVs_byinstitute.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =7, units = 'in')


# free y axis
p <- potu.c %>%
  # remove any nan or 0 read values
  filter(reads>0) %>%
  mutate(count = 1) %>%
  group_by(SampleID) %>%
  mutate(count = sum(count)) %>%
  ungroup() %>%
  distinct(SampleID, count) %>%
  right_join(meta_tab) %>%
  filter(site !='NA') %>%
  #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
  ggplot(aes(y=count, x=fct_relevel(Analyzing_Institute, institutes), color=fct_relevel(site, site_list))) +
  geom_boxplot() +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') +
  #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
  facet_wrap(fct_relevel(site, site_list)~.,  ncol=2, scales="free") +
  #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  labs(y='Number of ASVs', x='', color='site')
p

filename = paste(plot_dir, marker,'_totASVs_bysite.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =7, units = 'in')

# Plot total number of unique taxa ------------------------------------------

taxa_merged <- potu.c %>%
  filter(reads>0) %>%
  left_join(species_label %>% select(species_join, ASV)) %>%
  group_by(species_join, SampleID) %>%
  mutate(reads = sum(reads)) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  filter(reads>0) %>%
  distinct(species_join, SampleID, .keep_all = TRUE)
taxa_merged 

# by institute
p <- taxa_merged %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(count = sum(count)) %>%
  ungroup() %>%
  distinct(SampleID, count) %>%
  right_join(meta_tab) %>%
  filter(site !='NA') %>%
  #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
  ggplot(aes(y=count, x=fct_relevel(site, site_list), color=fct_relevel(site, site_list))) +
  geom_boxplot() +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') +
  #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
  facet_wrap(fct_relevel(Analyzing_Institute, institutes)~.,  ncol=2, scales="free") +
  #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  labs(y='Number of unique taxa', x='', color='site')+
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))
#theme(axis.text.x = element_text(angle=45, vjust= 0.3, hjust = 0.5))
p

filename = paste(plot_dir, marker,'_totUniqTaxa_byinstitute.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =7, units = 'in')

# by site
p <- taxa_merged %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(count = sum(count)) %>%
  ungroup() %>%
  distinct(SampleID, count) %>%
  right_join(meta_tab) %>%
  filter(site !='NA') %>%
  #ggplot(aes(y=reads, x=fct_relevel(Analyzing_Institute, institutes), color=site)) +
  #ggplot(aes(y=count, x=fct_relevel(site, site_list), color=fct_relevel(site, site_list))) +
  ggplot(aes(y=count, x=fct_relevel(Analyzing_Institute, institutes), color=fct_relevel(Analyzing_Institute, institutes))) +
  geom_boxplot() +
  geom_point() +
  scale_color_brewer(palette = 'Set1') +
  #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
  facet_wrap(fct_relevel(site, site_list)~.,  ncol=2, scales="free") +
  #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  labs(y='Number of unique taxa', x='', color='institute')+
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))
#theme(axis.text.x = element_text(angle=45, vjust= 0.3, hjust = 0.5))
p

filename = paste(plot_dir, marker,'_totUniqTaxa_bysite.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =7, units = 'in')

# Plot reads vs ASV/Taxa unique # ---------------------------------------------
# merged by taxonomy - get number of unique taxa
df1 <- taxa_merged %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(count = sum(count)) %>%
  ungroup() %>%
  distinct(SampleID, count) %>%
  rename(uniqTax = count)

# sum ASVs
df2 <- potu.c %>%
  # remove any nan or 0 read values
  filter(reads>0) %>%
  mutate(count = 1) %>%
  group_by(SampleID) %>%
  mutate(count = sum(count)) %>%
  ungroup() %>%
  distinct(SampleID, count) %>%
  rename(numASVs = count)

# sum reads
df3 <- potu.c %>%
  filter(reads>0) %>%
  group_by(SampleID) %>%
  mutate(reads = sum(reads, na.rm=TRUE)) %>%
  ungroup() %>%
  distinct(SampleID, reads) 
  
p <- df1 %>% full_join(df2) %>%
  full_join(df3) %>%
  right_join(meta_tab) %>%
  #filter(site !='NA') %>%
  ggplot(aes(y=reads, x=uniqTax, color=fct_relevel(site, site_list))) +
  #geom_boxplot() +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') +
  #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
  facet_wrap(fct_relevel(Analyzing_Institute, institutes)~.,  ncol=2, scales="free") +
  #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  labs(y='Total Reads', x='Unique Taxonomic Annotations', color='site')+
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))

p

filename = paste(plot_dir, marker,'_ReadsvTaxa.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =8, units = 'in')

p <- df1 %>% full_join(df2) %>%
  full_join(df3) %>%
  right_join(meta_tab) %>%
  #filter(site !='NA') %>%
  ggplot(aes(y=reads, x=numASVs, color=fct_relevel(site, site_list))) +
  #geom_boxplot() +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') +
  #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
  facet_wrap(fct_relevel(Analyzing_Institute, institutes)~.,  ncol=2, scales="free") +
  #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
  labs(y='Total Reads', x='Total ASVs', color='site')+
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))

p

filename = paste(plot_dir, marker,'_ReadsvASVs.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =8, units = 'in')

# # just look at AWI and the replicates:
# p <- df1 %>% full_join(df2) %>%
#   full_join(df3) %>%
#   left_join(meta_tab) %>%
#   filter(site !='NA') %>%
#   filter(Analyzing_Institute == 'AWI') %>%
#   # Framstrait_01_0008 Framstrait_01_0049
#   mutate(plate_string = str_match(SampleID, '.*_.*_(.*)')[,2]) %>%
#   ggplot(aes(y=reads, x=uniqTax, color=fct_relevel(site, site_list), shape=plate_string)) +
#   #geom_boxplot() +
#   geom_point() +
#   scale_color_brewer(palette = 'Dark2') +
#   #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
#   #facet_wrap(fct_relevel(Analyzing_Institute, institutes)~.,  ncol=2, scales="free") +
#   #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
#   labs(y='Total Reads', x='Unique Taxonomic Annotations', color='site')+
#   theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))
# 
# p
# 
# filename = paste(plot_dir, marker,'_ReadsvTaxa_AWI.png', sep='')
# #print('Plot of top 20 genus average by month:')
# print(filename)
# ggsave(filename,height = 6, width =6, units = 'in')
# 
# p <- df1 %>% full_join(df2) %>%
#   full_join(df3) %>%
#   left_join(meta_tab) %>%
#   filter(site !='NA') %>%
#   filter(Analyzing_Institute == 'AWI') %>%
#   # Framstrait_01_0008 Framstrait_01_0049
#   mutate(plate_string = str_match(SampleID, '.*_.*_(.*)')[,2]) %>%
#   ggplot(aes(y=reads, x=numASVs, color=fct_relevel(site, site_list), shape=plate_string)) +
#   #geom_boxplot() +
#   geom_point() +
#   scale_color_brewer(palette = 'Dark2') +
#   #facet_grid(fct_relevel(site, site_list)~fct_relevel(Analyzing_Institute, institutes), scales="free") +
#   facet_wrap(fct_relevel(Analyzing_Institute, institutes)~.,  ncol=2, scales="free") +
#   #facet_grid(fct_relevel(Analyzing_Institute, institutes)~fct_relevel(site, site_list), scales="free")
#   labs(y='Total Reads', x='Total ASVs', color='site')+
#   theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))
# 
# p
# 
# filename = paste(plot_dir, marker,'_ReadsvASVs_AWI.png', sep='')
# #print('Plot of top 20 genus average by month:')
# print(filename)
# ggsave(filename,height = 6, width =6, units = 'in')


# Unique taxa per taxonomic level -----------------------

# merged by taxonomy - get number of unique taxa per level

df<- potu.c %>%
  full_join(tax.c) %>%
  filter(reads>0) %>%
  group_by(SampleID) %>%
  mutate(domain_count = n_distinct(domain)) %>%
  mutate(supergroup_count = n_distinct(supergroup)) %>%
  mutate(division_count = n_distinct(division)) %>%
  mutate(subdivision_count = n_distinct(subdivision)) %>%
  mutate(class_count = n_distinct(class)) %>%
  mutate(order_count = n_distinct(order)) %>%
  mutate(family_count = n_distinct(family)) %>%
  mutate(genus_count = n_distinct(genus)) %>%
  mutate(species_count = n_distinct(species)) %>%
  ungroup() %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  select(SampleID,domain_count, supergroup_count, division_count, subdivision_count, class_count, order_count, family_count, genus_count, species_count) %>%
  pivot_longer(-SampleID, values_to = "count", names_to = "tax_level")

t_levels = c('domain_count', 'supergroup_count', 'division_count', 'subdivision_count', 'class_count', 'order_count', 'family_count', 'genus_count', 'species_count')
p <- df %>%
  right_join(meta_tab) %>%
  #filter(site %in% c('BLOOMMOCK', 'EVENMOCK', 'NOC')) %>%
  filter(site %in% c('Western Channel')) %>%
  ggplot(aes(x=fct_relevel(tax_level, t_levels), y=count, color=Analyzing_Institute)) +
  geom_point() +
  geom_line(aes(group=interaction(Analyzing_Institute,site, replicateID)), alpha=0.5)+
  labs(y='Number of Unique Annotations', x='taxonomic level') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))
p  
filename = paste(plot_dir, marker,'_unique_taxa_perlevel.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =8, units = 'in')


# Unique taxa per upper level group per sample ----------------------------------

divisions = c('Alveolata', 'Chlorophyta', 'Haptophyta', 'Cryptophyta', 'Rhizaria', 'Stramenopiles')

df<- potu.c %>%
  full_join(tax.c) %>%
  filter(reads>0) %>%
  filter(division=='Alveolata') %>%
  group_by(SampleID) %>%
  mutate(domain_count = n_distinct(domain)) %>%
  mutate(supergroup_count = n_distinct(supergroup)) %>%
  mutate(division_count = n_distinct(division)) %>%
  mutate(subdivision_count = n_distinct(subdivision)) %>%
  mutate(class_count = n_distinct(class)) %>%
  mutate(order_count = n_distinct(order)) %>%
  mutate(family_count = n_distinct(family)) %>%
  mutate(genus_count = n_distinct(genus)) %>%
  mutate(species_count = n_distinct(species)) %>%
  ungroup() %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  select(SampleID,domain_count, supergroup_count, division_count, subdivision_count, class_count, order_count, family_count, genus_count, species_count) %>%
  pivot_longer(-SampleID, values_to = "count", names_to = "tax_level")

p <- df %>%
  right_join(meta_tab) %>%
  #filter(site %in% c('BLOOMMOCK', 'EVENMOCK', 'NOC')) %>%
  filter(site %in% c('Western Channel')) %>%
  ggplot(aes(x=fct_relevel(tax_level, t_levels), y=count, color=Analyzing_Institute)) +
  geom_point() +
  geom_line(aes(group=interaction(Analyzing_Institute,site, replicateID)), alpha=0.5)+
  labs(y='Number of Unique Annotations', x='taxonomic level') +
  ggtitle('Western Channel Alveolata') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))


p 

#Stramenopiles
df<- potu.c %>%
  full_join(tax.c) %>%
  filter(reads>0) %>%
  filter(division=='Stramenopiles') %>%
  group_by(SampleID) %>%
  mutate(domain_count = n_distinct(domain)) %>%
  mutate(supergroup_count = n_distinct(supergroup)) %>%
  mutate(division_count = n_distinct(division)) %>%
  mutate(subdivision_count = n_distinct(subdivision)) %>%
  mutate(class_count = n_distinct(class)) %>%
  mutate(order_count = n_distinct(order)) %>%
  mutate(family_count = n_distinct(family)) %>%
  mutate(genus_count = n_distinct(genus)) %>%
  mutate(species_count = n_distinct(species)) %>%
  ungroup() %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  select(SampleID,domain_count, supergroup_count, division_count, subdivision_count, class_count, order_count, family_count, genus_count, species_count) %>%
  pivot_longer(-SampleID, values_to = "count", names_to = "tax_level")

p <- df %>%
  right_join(meta_tab) %>%
  #filter(site %in% c('BLOOMMOCK', 'EVENMOCK', 'NOC')) %>%
  filter(site %in% c('Western Channel')) %>%
  ggplot(aes(x=fct_relevel(tax_level, t_levels), y=count, color=Analyzing_Institute)) +
  geom_point() +
  geom_line(aes(group=interaction(Analyzing_Institute,site, replicateID)), alpha=0.5)+
  labs(y='Number of Unique Annotations', x='taxonomic level') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))


p 

#Chlorophyta
df<- potu.c %>%
  full_join(tax.c) %>%
  filter(reads>0) %>%
  filter(division=='Chlorophyta') %>%
  group_by(SampleID) %>%
  mutate(domain_count = n_distinct(domain)) %>%
  mutate(supergroup_count = n_distinct(supergroup)) %>%
  mutate(division_count = n_distinct(division)) %>%
  mutate(subdivision_count = n_distinct(subdivision)) %>%
  mutate(class_count = n_distinct(class)) %>%
  mutate(order_count = n_distinct(order)) %>%
  mutate(family_count = n_distinct(family)) %>%
  mutate(genus_count = n_distinct(genus)) %>%
  mutate(species_count = n_distinct(species)) %>%
  ungroup() %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  select(SampleID,domain_count, supergroup_count, division_count, subdivision_count, class_count, order_count, family_count, genus_count, species_count) %>%
  pivot_longer(-SampleID, values_to = "count", names_to = "tax_level")

p <- df %>%
  right_join(meta_tab) %>%
  #filter(site %in% c('BLOOMMOCK', 'EVENMOCK', 'NOC')) %>%
  filter(site %in% c('Western Channel')) %>%
  ggplot(aes(x=fct_relevel(tax_level, t_levels), y=count, color=Analyzing_Institute)) +
  geom_point() +
  geom_line(aes(group=interaction(Analyzing_Institute,site, replicateID)), alpha=0.5)+
  labs(y='Number of Unique Annotations', x='taxonomic level') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))


p 
#Haptophyta
df<- potu.c %>%
  full_join(tax.c) %>%
  filter(reads>0) %>%
  filter(division=='Haptophyta') %>%
  group_by(SampleID) %>%
  mutate(domain_count = n_distinct(domain)) %>%
  mutate(supergroup_count = n_distinct(supergroup)) %>%
  mutate(division_count = n_distinct(division)) %>%
  mutate(subdivision_count = n_distinct(subdivision)) %>%
  mutate(class_count = n_distinct(class)) %>%
  mutate(order_count = n_distinct(order)) %>%
  mutate(family_count = n_distinct(family)) %>%
  mutate(genus_count = n_distinct(genus)) %>%
  mutate(species_count = n_distinct(species)) %>%
  ungroup() %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  select(SampleID,domain_count, supergroup_count, division_count, subdivision_count, class_count, order_count, family_count, genus_count, species_count) %>%
  pivot_longer(-SampleID, values_to = "count", names_to = "tax_level")

p <- df %>%
  right_join(meta_tab) %>%
  #filter(site %in% c('BLOOMMOCK', 'EVENMOCK', 'NOC')) %>%
  filter(site %in% c('Western Channel')) %>%
  ggplot(aes(x=fct_relevel(tax_level, t_levels), y=count, color=Analyzing_Institute)) +
  geom_point() +
  geom_line(aes(group=interaction(Analyzing_Institute,site, replicateID)), alpha=0.5)+
  labs(y='Number of Unique Annotations', x='taxonomic level') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))


p 


# species level across divisions
#divisions = c('Alveolata', 'Chlorophyta', 'Haptophyta', 'Cryptophyta', 'Rhizaria', 'Stramenopiles')
divisions = c('Alveolata', 'Chlorophyta', 'Cryptophyta', 'Haptophyta', 'Rhizaria', 'Stramenopiles')
df<- potu.c %>%
  #full_join(tax.c) %>%
  left_join(species_label) %>%
  filter(reads>0) %>%
  select(SampleID, division, species_join) %>%
  group_by(SampleID, division) %>%
  mutate(sp_count = n_distinct(species_join)) %>%
  ungroup() %>%
  distinct(SampleID,division,sp_count) %>%
  filter(division %in% divisions)

p <- df %>%
  right_join(meta_tab) %>%
  #filter(site %in% c('BLOOMMOCK', 'EVENMOCK', 'NOC')) %>%
  filter(site %in% c('Western Channel')) %>%
  ggplot(aes(x=fct_relevel(division, divisions), y=sp_count, color=fct_relevel(Analyzing_Institute, institutes))) +
  geom_point() +
  geom_line(aes(group=interaction(Analyzing_Institute,site, replicateID)), alpha=0.5)+
  scale_color_brewer(palette = 'Set1') +
  labs(y='Number of Unique Annotations', x='taxonomic division') +
  ggtitle('Unrarefied unique annotations in Western Channel') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))


p 

filename = paste(plot_dir, marker,'_unique_taxa_westernchannel_bydivision.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =8, units = 'in')

p <- df %>%
  right_join(meta_tab) %>%
  filter(site %in% c('BLOOMMOCK', 'EVENMOCK')) %>%
  #filter(site %in% c('Western Channel')) %>%
  ggplot(aes(x=fct_relevel(division, divisions), y=sp_count, color=fct_relevel(Analyzing_Institute, institutes), fill =fct_relevel(Analyzing_Institute, institutes))) +
  geom_point() +
  #geom_bar(position='dodge',stat='identity') +
  #geom_bar() +
  #geom_boxplot() +
  geom_line(aes(group=interaction(Analyzing_Institute,site, replicateID)), alpha=0.5)+
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  labs(y='Number of Unique Annotations', x='taxonomic division') +
  ggtitle('Unrarefied unique annotations in controls') +
  theme(axis.text.x = element_text(angle=45, vjust= 1, hjust = 1))


p 
filename = paste(plot_dir, marker,'_unique_taxa_controls_bydivision.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 8, width =8, units = 'in')

# Make Phyloseq Object ----------------------------------------------------

library(phyloseq)

#ASV table
print('ASV table')
asv_path = paste(data_directory,'pr2_reassigned/phytoplankton_only/', prefix, '_asv_PR2_50_filtered.csv', sep='')
print(asv_path)

#taxa table
print('taxa table')
taxa_path = filepath = paste(data_directory, 'pr2_reassigned/phytoplankton_only/',prefix, '_taxa_PR2_50_filtered.csv', sep='')
print(taxa_path)

#metadata table
print('metadata table')
samp_path = paste(data_directory,'merged_original_taxonomy/', prefix, '_meta_merged.csv', sep='')
print(samp_path)

#Make phyloseq object
full_df <- merge_phyloseq(
  otu_table(as.matrix(read.csv(file = asv_path, row.names = 1,check.names = FALSE)), taxa_are_rows = TRUE),
  tax_table(as.matrix(read.csv(file = taxa_path, row.names = 1))),
  sample_data(read.csv(file = samp_path, row.names = 1)))
full_df
 


merged_df <- subset_samples(full_df, !(Collecting_Institute %in% c('BLOOMMOCK', 'EVENMOCK')) )
# remove duplicate AWI samples
merged_df <- subset_samples(merged_df, !grepl("0008", sample_names(merged_df)))
# remove blank
#blank_01_0049
merged_df <- subset_samples(merged_df, !grepl("blank_01_0049", sample_names(merged_df)))
NOC_only <- subset_samples(full_df, (Collecting_Institute %in% c('NOC')) )

## rarefy ---------------------------

print(min(sample_sums(merged_df)))
print(max(sample_sums(merged_df)))
print(mean(sample_sums(merged_df)))
print(merged_df)
# Remove samples with <500 reads (like DEICODE)
oBiom_r = prune_samples(sample_sums(merged_df)>=500, merged_df)
#set random seed for reproducibility later
print(min(sample_sums(oBiom_r)))
oBiom_r = rarefy_even_depth(merged_df, sample.size = min(sample_sums(oBiom_r)), rngseed = 678, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
print(oBiom_r)

# NOC only
NOC_oBiom_r = prune_samples(sample_sums(NOC_only)>=500, NOC_only)
#set random seed for reproducibility later
print(min(sample_sums(NOC_oBiom_r)))
NOC_oBiom_r = rarefy_even_depth(NOC_only, sample.size = min(sample_sums(NOC_oBiom_r)), rngseed = 678, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
print(NOC_oBiom_r)

## save rarefied data ------------------
#create matrices
rvegan_matrix = vegan_otu(oBiom_r)
rotu_tab = otu_table(oBiom_r)
rtax_tab = tax_table(oBiom_r)
# rsamp_dat = sample_data(oBiom_r)


#save rarefied data
#save to csv file
file = paste(data_directory,'/pr2_reassigned/phytoplankton_only/rarefied_data/',marker,'_Rarefied_asv_table.csv', sep="")
write.csv(rotu_tab, file)
file = paste(data_directory,'/pr2_reassigned/phytoplankton_only/rarefied_data/',marker,'_Rarefied_taxa_table.csv', sep="")
write.csv(rtax_tab, file)
# file = paste(data_directory,'/pr2_reassigned/phytoplankton_only/rarefied_data/',marker,'_Rarefied_meta_table.csv', sep="")
# write.csv(rsamp_dat_tibble, file)


# diversity indices --------------------
library(ggplot2)
p_oBiom <- plot_richness(oBiom_r,x="site", measures=c("Shannon", "Chao1", "Simpson"), nrow=1, color='Analyzing_Institute', shape='Analyzing_Institute') +
  theme(text = element_text(size=16), strip.text.x = element_text(size = 16), axis.text = element_text(size = rel(1), colour = "black")) +
  xlab("Site") + 
  ggtitle(paste("Community alpha diversity: ",marker,sep=""))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_brewer(palette = "Spectral")+scale_color_brewer(palette = "Dark2")
p_oBiom 
filename = paste(plot_dir, marker,'_diversity_indices_phyloseq.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

p_oBiom <- plot_richness(oBiom_r,x="site", measures=c("Chao1"), nrow=1, color='Analyzing_Institute', shape='Analyzing_Institute') +
  geom_boxplot() +
  theme(text = element_text(size=16), strip.text.x = element_text(size = 16), axis.text = element_text(size = rel(1), colour = "black")) +
  xlab("Site") + 
  ggtitle(paste("Community alpha diversity: ",marker,sep=""))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_brewer(palette = "Spectral")+scale_color_brewer(palette = "Dark2")
p_oBiom 
filename = paste(plot_dir, marker,'_diversity_Chao1_phyloseq.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')



# NOC only
p_oBiom <- plot_richness(NOC_oBiom_r,x="Analyzing_Institute", measures=c("Shannon", "Chao1", "Simpson"), nrow=1, color='Analyzing_Institute', shape='Analyzing_Institute') +
  theme(text = element_text(size=16), strip.text.x = element_text(size = 16), axis.text = element_text(size = rel(1), colour = "black")) +
  #xlab("Site") + 
  ggtitle(paste("Community alpha diversity: ",marker,sep=""))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_brewer(palette = "Spectral")+scale_color_brewer(palette = "Dark2")
p_oBiom 
filename = paste(plot_dir, marker,'_diversity_indices_phyloseq_NOC.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')



###  ranacapa rarefaction curves ----------------
library(ranacapa)
# p <- ggrare(merged_df, step = 1000, label = "Analyzing_Institute", color = "Collecting_Institute",
#             plot = TRUE, parallel = FALSE, se = FALSE)
# p


p <- ggrare(merged_df, step = 1000, label = "Collecting_Institute", color = "Analyzing_Institute",
            plot = TRUE, parallel = FALSE, se = FALSE)
p


filename = paste(plot_dir, marker,'_diversity_rarefaction_curves.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 5, width =12, units = 'in')

### merge by taxonomy ----------------

#### family 
ps_taxglom <- tax_glom(physeq = merged_df, taxrank = "family")

nmds_ordination <- ordinate(ps_taxglom, method = "NMDS", distance = "bray")
p <- plot_ordination(ps_taxglom, nmds_ordination, color = "site", shape = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("Family-level NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_collecting_family.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

p <- plot_ordination(ps_taxglom, nmds_ordination, shape = "site", color = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("Family-level NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_analyzing_family.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

# Calculate distance matrix (e.g., Bray-Curtis)
dist_matrix <- phyloseq::distance(ps_taxglom, method = "bray")
metadata_df <- data.frame(sample_data(ps_taxglom))
# Perform PERMANOVA using adonis2 from vegan
# 'GroupVariable' is a categorical variable in your sample data
permanova_result <- vegan::adonis2(dist_matrix ~ site + Analyzing_Institute, data = metadata_df, by = "terms")
#permanova_result <- vegan::adonis2(dist_matrix ~ site, data = metadata_df, permutations = 99, strata = metadata_df$Analyzing_Institute)


# View the results
print(permanova_result)

##### Jaccard_family

# ps_taxglom <- tax_glom(physeq = merged_df, taxrank = "family")
# jsd, gower, w, canberra, kulczynski
nmds_ordination <- ordinate(ps_taxglom, method = "NMDS", distance = "jsd")
p <- plot_ordination(ps_taxglom, nmds_ordination, color = "site", shape = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("Family-level NMDS Ordination with jsd") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_collecting_family_jsd.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

p <- plot_ordination(ps_taxglom, nmds_ordination, shape = "site", color = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("Family-level NMDS Ordination with jsd") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_analyzing_family_jsd.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

# Calculate distance matrix (e.g., Bray-Curtis)
dist_matrix <- phyloseq::distance(ps_taxglom, method = "jsd")
metadata_df <- data.frame(sample_data(ps_taxglom))
# Perform PERMANOVA using adonis2 from vegan
# 'GroupVariable' is a categorical variable in your sample data
permanova_result <- vegan::adonis2(dist_matrix ~ site + Analyzing_Institute, data = metadata_df, by = "terms")
#permanova_result <- vegan::adonis2(dist_matrix ~ site, data = metadata_df, permutations = 99, strata = metadata_df$Analyzing_Institute)


# View the results
print(permanova_result)


#### order
ps_taxglom <- tax_glom(physeq = merged_df, taxrank = "order")

nmds_ordination <- ordinate(ps_taxglom, method = "NMDS", distance = "bray")
p <- plot_ordination(ps_taxglom, nmds_ordination, color = "site", shape = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_collecting_order.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

p <- plot_ordination(ps_taxglom, nmds_ordination, shape = "site", color = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_analyzing_order.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

#### genus
ps_taxglom <- tax_glom(physeq = merged_df, taxrank = "genus")

nmds_ordination <- ordinate(ps_taxglom, method = "NMDS", distance = "bray")
p <- plot_ordination(ps_taxglom, nmds_ordination, color = "site", shape = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_collecting_genus.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

p <- plot_ordination(ps_taxglom, nmds_ordination, shape = "site", color = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_analyzing_genus.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

nmds_ordination <- ordinate(ps_taxglom, method = "NMDS", distance = "jaccard")
p <- plot_ordination(ps_taxglom, nmds_ordination, color = "site", shape = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_collecting_genus_jaccard.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

p <- plot_ordination(ps_taxglom, nmds_ordination, shape = "site", color = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_analyzing_genus_jaccard.png', sep='')
#print('Plot of top 20 genus average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

#### species
ps_taxglom <- tax_glom(physeq = merged_df, taxrank = "species")

nmds_ordination <- ordinate(ps_taxglom, method = "NMDS", distance = "bray")
p <- plot_ordination(ps_taxglom, nmds_ordination, color = "site", shape = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_collecting_species.png', sep='')
#print('Plot of top 20 species average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

p <- plot_ordination(ps_taxglom, nmds_ordination, shape = "site", color = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_analyzing_species.png', sep='')
#print('Plot of top 20 species average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

nmds_ordination <- ordinate(ps_taxglom, method = "NMDS", distance = "jaccard")
p <- plot_ordination(ps_taxglom, nmds_ordination, color = "site", shape = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_collecting_species_jaccard.png', sep='')
#print('Plot of top 20 species average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

p <- plot_ordination(ps_taxglom, nmds_ordination, shape = "site", color = "Analyzing_Institute") +
  geom_point(size = 3) + # Adjust point size for better visualization
  stat_ellipse() + # Add ellipses to group samples by categories
  ggtitle("NMDS Ordination with Bray-Curtis Dissimilarity") +
  theme_bw() # Use a clean black and white theme

p
filename = paste(plot_dir, marker,'_nmds_color_analyzing_species_jaccard.png', sep='')
#print('Plot of top 20 species average by month:')
print(filename)
ggsave(filename,height = 6.5, width =8, units = 'in')

## upset plot by group ---------------
library(UpSetR)

tibb_rotu <-tibble::rownames_to_column(as.data.frame(rotu_tab), var = "ASV")
tibb_rotu <- tibble::as_tibble(tibb_rotu)
# convert to long format:
tibb_rotu %<>%
  tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
  group_by(SampleID) %>%
  mutate(per_tot = reads / sum(reads, na.rm=TRUE) *100) %>%
  ungroup() %>%
  arrange(-reads)

# environmental data - rarefied
env_data <- inner_join(meta_tab, tibb_rotu, by = c("SampleID")) %>% 
  left_join(species_label)

env_data 


#plot by plate
p1 <- tibb_rotu %>%
  left_join(meta_tab %>% select(SampleID, Analyzing_Institute, Collecting_Institute, site)) %>%
  left_join(species_label) %>%
  filter(Collecting_Institute == 'NOC') %>%
  # filter(depth <=300) %>%  # try to remove depth effects
  # filter(reads >0) %>%
  # filter(Phylum == 'Chordata') %>%
  # filter(Class == 'Actinopteri') %>%
  group_by(Analyzing_Institute, species_join) %>%
  mutate(reads=sum(reads)) %>%
  ungroup() %>%
  distinct(Analyzing_Institute, reads, species_join) %>%
  # Put a read limit in by taxa:
  group_by(Analyzing_Institute,species_join) %>%
  mutate(Total = sum(reads)) %>%
  ungroup() %>%
  filter(Total >0) %>%
  mutate(count = 1) %>%
  #filter(reads >50) %>%
  #unite(taxa,Kingdom, Phylum, Class, Order, Family, Genus, 
  #     Species) %>%
  #group_by(PlateID, taxa) %>%
  #count(taxa) %>%
  pivot_wider(
    #id_cols = taxa,
    id_cols = c(species_join),
    names_from = Analyzing_Institute,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  data.frame()

upset(p1)

#plot by plate - genus
p1 <- tibb_rotu %>%
  left_join(meta_tab %>% select(SampleID, Analyzing_Institute, Collecting_Institute, site)) %>%
  left_join(species_label) %>%
  filter(Collecting_Institute == 'NOC') %>%
  group_by(Analyzing_Institute, genus_join) %>%
  mutate(reads=sum(reads)) %>%
  ungroup() %>%
  distinct(Analyzing_Institute, reads, genus_join) %>%
  # Put a read limit in by taxa:
  group_by(Analyzing_Institute,genus_join) %>%
  mutate(Total = sum(reads)) %>%
  ungroup() %>%
  filter(Total >0) %>%
  mutate(count = 1) %>%
  pivot_wider(
    id_cols = c(genus_join),
    names_from = Analyzing_Institute,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  data.frame()

upset(p1)

# upset - family - complex --------
library(ComplexUpset)

p1 <- tibb_rotu %>%
  left_join(meta_tab %>% select(SampleID, Analyzing_Institute, Collecting_Institute, site)) %>%
  left_join(species_label) %>%
  filter(Collecting_Institute == 'NOC') %>%
  group_by(Analyzing_Institute, division, family_join) %>%
  mutate(reads=sum(reads)) %>%
  ungroup() %>%
  distinct(Analyzing_Institute, division, reads, family_join) %>%
  # # Put a read limit in by taxa:
  # group_by(Analyzing_Institute,family_join) %>%
  # mutate(Total = sum(reads)) %>%
  # ungroup() %>%
  filter(reads >0) %>%
  mutate(count = 1) %>%
  pivot_wider(
    id_cols = c(family_join, division),
    names_from = Analyzing_Institute,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  data.frame()

insts <- c('SBR', 'MBARI', 'UDAL', 'AOML', 'AWI')
# p <- ComplexUpset::upset(p1, intersect=insts)

p <- ComplexUpset::upset(
  p1,
  insts,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=division)
    ) + scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
  ),
  width_ratio=0.1
)


filename = paste(plot_dir, marker,'_upset_family.png', sep='')
#print('Plot of top 20 species average by month:')
print(filename)
ggsave(filename,height = 8, width =12, units = 'in')

# upset - genus - complex --------
# library(ComplexUpset)

p1 <- tibb_rotu %>%
  left_join(meta_tab %>% select(SampleID, Analyzing_Institute, Collecting_Institute, site)) %>%
  left_join(species_label) %>%
  filter(Collecting_Institute == 'NOC') %>%
  group_by(Analyzing_Institute, division, family_join) %>%
  mutate(reads=sum(reads)) %>%
  ungroup() %>%
  distinct(Analyzing_Institute, division, reads, family_join) %>%
  # # Put a read limit in by taxa:
  # group_by(Analyzing_Institute,family_join) %>%
  # mutate(Total = sum(reads)) %>%
  # ungroup() %>%
  filter(reads >0) %>%
  mutate(count = 1) %>%
  pivot_wider(
    id_cols = c(family_join, division),
    names_from = Analyzing_Institute,
    values_from = count,
    values_fill = list(count = 0)
  ) %>%
  data.frame()

insts <- c('SBR', 'MBARI', 'UDAL', 'AOML', 'AWI')
# p <- ComplexUpset::upset(p1, intersect=insts)

p <- ComplexUpset::upset(
  p1,
  insts,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=division)
    ) + scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
  ),
  width_ratio=0.1
)


filename = paste(plot_dir, marker,'_upset_family.png', sep='')
#print('Plot of top 20 species average by month:')
print(filename)
ggsave(filename,height = 8, width =14, units = 'in')


