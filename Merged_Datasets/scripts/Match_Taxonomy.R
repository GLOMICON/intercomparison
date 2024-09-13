# Apply Taxize to standardize taxonomy between datasets

# Set Directories ----------------------

#locations to store files:
#data files
data_directory <- './Merged_Datasets/data/'
#figures
plot_dir <- './Merged_Datasets/figures/taxonomy/'

# Load Libraries ----------------------
library(tidyverse)
library(readr)
library(taxize)
library(ggthemes)
library(forcats)
# library(ggbump)

# Functions --------------------------

make_compositional <- function(df) {
  df %<>%
    tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
    group_by(SampleID) %>%
    mutate(per_tot = reads / sum(reads,na.rm=TRUE ) *100) %>%
    ungroup() %>%
    arrange(-reads)
  return(df)
}



# Load ASV to Institution Table ------------------
filename <- paste(data_directory, 'GLOMICON_seq_merged.csv', sep='')

seq_tab <- read_csv(filename)

filename <- paste(data_directory, 'GLOMICON_asv_merged.csv', sep='')

asv_tab <- read_csv(filename)

filename <- paste(data_directory, 'GLOMICON_meta_merged.csv', sep='')

meta_tab <- read_csv(filename)

potu_tab <- make_compositional(asv_tab)

# Load Taxonomy Table -----------------
filename <- paste(data_directory, 'GLOMICON_taxa_merged.csv', sep='')

taxa_tab <- read_csv(filename)

# MBARI ASVs already have NCBI taxonomy -----------------

# # run them anyway - comment out
# MBARI_tax <- seq_tab %>% full_join(taxa_tab) %>%
#   filter(Analyzing_Institute=='MBARI') 
# 
# other_tax <- seq_tab %>% full_join(taxa_tab) %>%
#   filter(Analyzing_Institute!='MBARI')


# Look at Species Level Taxonomy  ---------------------
# taxize::use_entrez()


GLOMICON_species_key <- seq_tab %>%
  full_join(taxa_tab) %>%
  mutate(Species_edited = Species) %>%
  mutate(Species_edited = str_replace(Species_edited, 'Maxillopoda', 'Hexanauplia')) %>%
  #	Gonyaulax_spinifera - make sure this name survives editing to remove _sp.
  # R you have to double escape regex characters - \\
  mutate(Species_edited = str_replace(Species_edited, '_X*_sp\\.', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, '_X*_', '')) %>%
  #mutate(Species_edited = str_replace(Species_edited, '_XX_sp\\.', '')) %>%
  #mutate(Species_edited = str_replace(Species_edited, '_XXX_sp\\.', '')) %>%
  #mutate(Species_edited = str_replace(Species_edited, '_XXXX_sp\\.', '')) %>%
  #mutate(Species_edited = str_replace(Species_edited, '_XXXXXX_sp\\.', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, 'XXXXX', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, 'XXXXXX', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, ':nucl', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, '_sp\\.', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, '_sp1', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, '_Group-1', '')) %>%
  mutate(Species_edited = str_replace(Species_edited,'_Clade.*', '')) %>%
  mutate(Species_edited = str_replace(Species_edited,'_Group.*', '')) %>%
  mutate(Species_edited = str_replace(Species_edited,'-Group.*', '')) %>%
  mutate(Species_edited = str_replace(Species_edited,'-lineage', '')) %>%
  mutate(Species_edited = str_replace(Species_edited,'_[1-9][a-z]*', '')) %>%
  mutate(Species_edited = str_replace(Species_edited,'_[a-z][1-9]', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, 'Basal', '')) %>%
  mutate(Species_edited = str_replace(Species_edited, '_', ' ')) %>%
  arrange(Species_edited)

GLOMICON_species <-GLOMICON_species_key %>%
  distinct(Species_edited, .keep_all = FALSE)

#batch
sp1 <- GLOMICON_species$Species_edited[1:600]
sp2 <- GLOMICON_species$Species_edited[601:900]
sp3 <- GLOMICON_species$Species_edited[901:1200]
sp4 <- GLOMICON_species$Species_edited[1201:1640]

#look up NCBI names

# passing in a long list of names causes error on NCBI side, batch instead:

# taxize_options(ncbi_sleep = 0.5) # just as an e.g., you can adjust this time
taxize_options(ncbi_sleep = 0.5)
out1 <- classification(sp1, db = 'ncbi', batch_size=5)
out2 <- classification(sp2, db = 'ncbi', batch_size=5)
out3 <- classification(sp3, db = 'ncbi', batch_size=5)
out4 <- classification(sp4, db = 'ncbi', batch_size=5)


# out1
out_test <- out1[!is.na(out1)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab1 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  filter(rank %in% c("kingdom","phylum","class","order","family","genus", "species")) %>% 
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class, order, family, genus, species)

# out2
out_test <- out2[!is.na(out2)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab2 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  filter(rank %in% c("kingdom","phylum","class","order","family","genus", "species")) %>% 
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class, order, family, genus, species)


# out3
out_test <- out3[!is.na(out3)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab3 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  filter(rank %in% c("kingdom", "phylum","class","order","family","genus", "species")) %>% 
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class, order, family, genus, species)

# out4
out_test <- out4[!is.na(out4)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab4 <- tibble(names = names(out_test), out_test) %>%
  unnest(cols = c(out_test)) %>%
  filter(rank %in% c("kingdom","phylum","class","order","family","genus", "species")) %>%
  select(-id) %>%
  spread(rank, name) %>%
  select(sci_name = names, kingdom, phylum, class, order, family, genus, species)

#merge taxa tibbles - only has species names with matches.
df <- full_join(test_tax_tab1,test_tax_tab2) %>%
  full_join(test_tax_tab3) %>%
  full_join(test_tax_tab4)

# original names are in column 'sci_name' in df
# In GLOMICON_species_key, edited species column name is in 'Species_edited'

# merge to original dataset - pull out things without matches.
# matching by 'Species_edited' name; became 'sci_name' in output of taxize
# create 'found' column to easily separate out which ASVs have results
test_sp <- full_join(df%>%mutate(found=1),GLOMICON_species_key %>% mutate(sci_name = Species_edited))
test_unassigned_sp <- test_sp %>%
  filter(is.na(found)==TRUE)

# Look at Genus Level -----------------------------------------

# Look for Genera matches among missing taxa:
GLOMICON_genera_key <- test_unassigned_sp %>%
  mutate(Genus_edited = Genus) %>%
  mutate(Genus_edited = str_replace(Genus_edited, 'not_provided', '')) %>%
  mutate(Genus_edited = str_replace(Genus_edited, 'Basal_', '')) %>%
  mutate(Genus_edited = str_replace(Genus_edited, '_.*', '')) %>%
  mutate(Genus_edited = str_replace(Genus_edited, '-.*', ''))

GLOMICON_genera <- GLOMICON_genera_key %>%
  distinct(Genus_edited)

gen1 <- GLOMICON_genera$Genus_edited[1:428]


out1_gen1 <- classification(gen1, db = 'ncbi', batch_size=5)


# out1
out_test <- out1_gen1[!is.na(out1_gen1)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab_gen1 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  filter(rank %in% c("kingdom","phylum","class","order","family","genus", "species")) %>% 
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class, order, family, genus)

# join with previous information
test_gen <- full_join(test_tax_tab_gen1%>%mutate(found=1),GLOMICON_genera_key %>% mutate(sci_name = Genus_edited) %>%
                        select(-kingdom, -phylum, -class, -order, -family, -genus, -species, -found), by='sci_name')
test_gen_unassigned <- test_gen %>%
  filter(is.na(found)==TRUE) #%>%

# Look at Family Level --------------------------

# Look for Family matches among missing taxa:
GLOMICON_family_key <- test_gen_unassigned %>%
  mutate(Family_edited = Family) %>%
  mutate(Family_edited = str_replace(Family_edited, 'not_provided', '')) %>%
  mutate(Family_edited = str_replace(Family_edited, 'Basal_', '')) %>%
  mutate(Family_edited = str_replace(Family_edited, '_.*', '')) %>%
  mutate(Family_edited = str_replace(Family_edited, '-.*', ''))

GLOMICON_family <- GLOMICON_family_key %>%
  distinct(Family_edited)

fam1 <- GLOMICON_family$Family_edited[1:301]


out1_fam1 <- classification(fam1, db = 'ncbi', batch_size=5)


# out1
out_test <- out1_fam1[!is.na(out1_fam1)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab_fam1 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  filter(rank %in% c("kingdom","phylum","class","order","family","genus", "species")) %>% 
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class, order, family)

# join with previous information
test_fam <- full_join(test_tax_tab_fam1%>%mutate(found=1),GLOMICON_family_key %>% mutate(sci_name = Family_edited) %>%
                        select(-kingdom, -phylum, -class, -order, -family, -genus,  -found), by='sci_name')
test_fam_unassigned <- test_fam %>%
  filter(is.na(found)==TRUE) 

# Look at Order Level --------------------------

# Look for Order matches among missing taxa:
GLOMICON_order_key <- test_fam_unassigned %>%
  mutate(Order_edited = Order) %>%
  mutate(Order_edited = str_replace(Order_edited, 'not_provided', '')) %>%
  mutate(Order_edited = str_replace(Order_edited, 'Basal_', '')) %>%
  #Sagenista clade getting misID'd as wasp genera
  mutate(Order_edited = str_replace(Order_edited, 'Sagenista', 'Bigyra')) %>%
  mutate(Order_edited = str_replace(Order_edited, '_.*', '')) %>%
  mutate(Order_edited = str_replace(Order_edited, '-.*', ''))

GLOMICON_order <- GLOMICON_order_key %>%
  distinct(Order_edited)

ord1 <- GLOMICON_order$Order_edited[1:133]


out1_ord1 <- classification(ord1, db = 'ncbi', batch_size=5)


# out1
out_test <- out1_ord1[!is.na(out1_ord1)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab_ord1 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  filter(rank %in% c("kingdom","phylum","class","order","family","genus", "species")) %>% 
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class, order)

# join with previous information
test_ord <- full_join(test_tax_tab_ord1%>%mutate(found=1),GLOMICON_order_key %>% mutate(sci_name = Order_edited) %>%
                        select(-kingdom, -phylum, -class, -order, -family, -found), by='sci_name')
test_ord_unassigned <- test_ord %>%
  filter(is.na(found)==TRUE) 

# protist group 'Sagenista' getting misID as wasp genus Sagenista

# Look at Class Level --------------------------

# Look for Class matches among missing taxa:
GLOMICON_class_key <- test_ord_unassigned %>%
  mutate(Class_edited = Class) %>%
  mutate(Class_edited = str_replace(Class_edited, 'not_provided', '')) %>%
  mutate(Class_edited = str_replace(Class_edited, 'Basal_', '')) %>%
  #protist Sagenista clade getting misID'd as wasp genera - bump up to Phylum
  #mutate(Class_edited = str_replace(Class_edited, 'Sagenista', 'Bigyra')) %>%
  mutate(Class_edited = str_replace(Class_edited, '_.*', '')) %>%
  mutate(Class_edited = str_replace(Class_edited, '-.*', ''))

GLOMICON_class <- GLOMICON_class_key %>%
  distinct(Class_edited)

cla1 <- GLOMICON_class$Class_edited[1:76]


out1_cla1 <- classification(cla1, db = 'ncbi', batch_size=5)


# out1
out_test <- out1_cla1[!is.na(out1_cla1)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab_cla1 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  filter(rank %in% c("kingdom","phylum","class","order","family","genus", "species")) %>% 
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class)

# join with previous information
test_cla <- full_join(test_tax_tab_cla1%>%mutate(found=1),GLOMICON_class_key %>% mutate(sci_name = Class_edited) %>%
                        select(-kingdom,-phylum, -class, -order,  -found), by='sci_name')
test_cla_unassigned <- test_cla %>%
  filter(is.na(found)==TRUE) 

# Look at Phylum Level --------------------------

# starting to run into problems with higher-level annotations.

# Look for Phylum matches among missing taxa:
GLOMICON_phylum_key <- test_ord_unassigned %>%
  mutate(Phylum_edited = Phylum) %>%
  mutate(Phylum_edited = str_replace(Phylum_edited, 'not_provided', '')) %>%
  mutate(Phylum_edited = str_replace(Phylum_edited, 'Basal_', '')) %>%
  #protist Sagenista clade getting misID'd as wasp genera - bump up to Phylum
  #mutate(Phylum_edited = str_replace(Phylum_edited, 'Sagenista', 'Bigyra')) %>%
  mutate(Phylum_edited = str_replace(Phylum_edited, '_.*', '')) %>%
  mutate(Phylum_edited = str_replace(Phylum_edited, '-.*', ''))

GLOMICON_phylum <- GLOMICON_phylum_key %>%
  distinct(Phylum_edited)

phy1 <- GLOMICON_phylum$Phylum_edited[1:52]


out1_phy1 <- classification(phy1, db = 'ncbi', batch_size=5)


# out1
out_test <- out1_phy1[!is.na(out1_phy1)]
tr <- class2tree(out_test)
plot(tr, no.margin = TRUE)

test_tax_tab_phy1 <- tibble(names = names(out_test), out_test) %>% 
  unnest(cols = c(out_test)) %>% 
  #filter(rank %in% c("superkingdom","kingdom","clade","phylum","class","order","family","genus", "species")) %>% 
  filter(rank %in% c("kingdom","phylum","class")) %>%
  select(-id) %>% 
  spread(rank, name) %>% 
  select(sci_name = names, kingdom, phylum, class)

# join with previous information
test_phy <- full_join(test_tax_tab_phy1%>%mutate(found=1),GLOMICON_phylum_key %>% mutate(sci_name = Phylum_edited) %>%
                        select(-kingdom,-phylum, -class, -order,  -found), by='sci_name')
test_phy_unassigned <- test_cla %>%
  filter(is.na(found)==TRUE) 


# Merge Information Together -------------------------
#distinct keeps first row (sort so first row is found=1, if available)
merged_df <- test_sp %>%
  full_join(test_gen) %>%
  full_join(test_fam) %>%
  full_join(test_ord) %>%
  full_join(test_cla) %>%
  full_join(test_phy) %>%
  arrange(ASV, found) %>%
  distinct(ASV, .keep_all = TRUE)
  
# export
filename <- paste(data_directory, 'GLOMICON_taxa_merged_updated_raw.csv', sep='')
#filename <- '/Users/kpitz/github/GLOMICON/intercomparison/Merged_Datasets/data/taxize_taxonomy/GLOMICON_taxa_merged_updated.csv'
print(filename)
write.csv(merged_df, filename)

# Get into right format for taxa table moving forward:

taxa_updated_df <- merged_df %>%
  select(ASV,kingdom, phylum, class, order, family, genus, species) %>%
  rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family= family, Genus = genus, Species = species)

# Format asv sequence table and taxa table; export -----------------

# # need to join back with MBARI data
# 
# updated_taxa_tab <- merged_df %>%
#   #select(-sci_name) %>%
#   select(ASV,kingdom, phylum, class, order, family, genus, species) %>%
#   rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family= family, Genus = genus, Species = species) %>%
#   bind_rows(MBARI_tax %>% select(-sequence, -Analyzing_Institute, -Domain, -taxonomy, -Supergroup))
#   #full_join(MBARI_tax, by='ASV')
# 
# 

#just check it matches with original seq_tab
test <-taxa_updated_df  %>% mutate(new_df = 1)%>%full_join(seq_tab, by='ASV')

# export
filename <- paste(data_directory, 'GLOMICON_taxa_merged_updated.csv', sep='')
#filename <- '/Users/kpitz/github/GLOMICON/intercomparison/Merged_Datasets/data/taxize_taxonomy/GLOMICON_taxa_merged_updated.csv'
print(filename)
write.csv(taxa_updated_df , filename)

# # Prelim look at data --------------
# # combine with seq_tab:
# 
# tot_reads = asv_tab %>%
#   tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
#   group_by(ASV) %>%
#   mutate(total_reads = sum(reads,na.rm=TRUE )) %>%
#   ungroup() %>%
#   distinct(ASV, .keep_all = TRUE) %>%
#   select(ASV, total_reads)
# 
# test <- asv_tab %>%
#   tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
#   full_join(seq_tab %>%select(ASV, Analyzing_Institute)) %>%
#   mutate(count=1) %>%
#   group_by(Analyzing_Institute) %>%
#   mutate(total_reads = sum(reads,na.rm=TRUE )) %>%
#   mutate(total_ASVs = sum(count, na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, .keep_all = TRUE) %>%
#   select(Analyzing_Institute, total_reads, total_ASVs) %>%
#   mutate(ASVs_per_read = total_ASVs/total_reads)
# 
# # look at NOC or UDalhousie where have equal replicates across all three analyzing institutes
# test <- potu_tab %>%
#   full_join(meta_tab %>% mutate(SampleID = sample_name)) %>%
#   filter(Collecting_Institute == 'NOC') %>%
#   group_by(Analyzing_Institute, ASV) %>%
#   mutate(sum_pertot = sum(per_tot)) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, ASV, .keep_all = TRUE) %>%
#   filter(sum_pertot >0) %>%
#   arrange(-sum_pertot) %>%
#   select(ASV, Analyzing_Institute, sum_pertot ) %>%
#   full_join(updated_taxa_tab %>% select(ASV, Phylum, Class, Order, Family, Genus, Species))
# 
# # number of ASVs per class:
# test <- potu_tab %>%
#   full_join(meta_tab %>% mutate(SampleID = sample_name)) %>%
#   filter(Collecting_Institute == 'NOC') %>%
#   group_by(Analyzing_Institute, ASV) %>%
#   mutate(sum_pertot = sum(per_tot)) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, ASV, .keep_all = TRUE) %>%
#   filter(sum_pertot >0) %>%
#   arrange(-sum_pertot) %>%
#   select(ASV, Analyzing_Institute, sum_pertot ) %>%
#   mutate(count=1) %>%
#   full_join(updated_taxa_tab %>% select(ASV, Phylum, Class, Order, Family, Genus, Species)) %>%
#   group_by(Analyzing_Institute, Class) %>%
#   mutate(sum_count = sum(count,na.rm=TRUE )) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, Class, .keep_all = TRUE) %>%
#   select(Analyzing_Institute, Class,sum_count) %>%
#   arrange(Class,Analyzing_Institute)
# 
# # number of ASVs per Phylum:
# test <- potu_tab %>%
#   full_join(meta_tab %>% mutate(SampleID = sample_name)) %>%
#   filter(Collecting_Institute %in% c('NOC', 'UDalhousie')) %>%
#   group_by(Analyzing_Institute, ASV) %>%
#   mutate(sum_pertot = sum(per_tot)) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, ASV, .keep_all = TRUE) %>%
#   filter(sum_pertot >0) %>%
#   arrange(-sum_pertot) %>%
#   select(ASV, Analyzing_Institute, sum_pertot ) %>%
#   mutate(count=1) %>%
#   full_join(updated_taxa_tab %>% select(ASV, Phylum, Class, Order, Family, Genus, Species)) %>%
#   group_by(Analyzing_Institute, Phylum) %>%
#   mutate(sum_count = sum(count,na.rm=TRUE )) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, Phylum, .keep_all = TRUE) %>%
#   select(Analyzing_Institute, Phylum,sum_count) %>%
#   arrange(Phylum,Analyzing_Institute)
# 
# 
# 
# 
# # look at NOC or UDalhousie where have equal replicates across all three analyzing institutes
# # total reads/ total asvs
# test <- potu_tab %>%
#   full_join(meta_tab %>% mutate(SampleID = sample_name)) %>%
#   filter(Collecting_Institute == 'NOC') %>%
#   mutate(count=1) %>%
#   filter(per_tot >0) %>%
#   group_by(Analyzing_Institute) %>%
#   mutate(total_reads = sum(reads,na.rm=TRUE )) %>%
#   mutate(total_ASVs = sum(count, na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, .keep_all = TRUE) %>%
#   select(Analyzing_Institute, total_reads, total_ASVs) %>%
#   mutate(ASVs_per_read = total_ASVs/total_reads)
# 
# test <- potu_tab %>%
#   full_join(meta_tab %>% mutate(SampleID = sample_name)) %>%
#   filter(Collecting_Institute == 'UDalhousie') %>%
#   mutate(count=1) %>%
#   filter(per_tot >0) %>%
#   group_by(Analyzing_Institute) %>%
#   mutate(total_reads = sum(reads,na.rm=TRUE )) %>%
#   mutate(total_ASVs = sum(count, na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(Analyzing_Institute, .keep_all = TRUE) %>%
#   select(Analyzing_Institute, total_reads, total_ASVs) %>%
#   mutate(ASVs_per_read = total_ASVs/total_reads)
# 
# 
# test <- seq_tab %>%
#   full_join(test_sp) %>%
#   select(class, Analyzing_Institute, ASV) %>%
#   mutate(count=1) %>%
#   full_join(tot_reads) %>%
#   group_by(class, Analyzing_Institute) %>%
#   mutate(sum_count = sum(count,na.rm=TRUE)) %>%
#   mutate(sum_reads = sum(total_reads,na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(class, Analyzing_Institute, .keep_all = TRUE) %>%
#   arrange(class, Analyzing_Institute)
# 
# test <- seq_tab %>%
#   full_join(test_sp) %>%
#   select(phylum, Analyzing_Institute, ASV) %>%
#   mutate(count=1) %>%
#   full_join(tot_reads) %>%
#   group_by(phylum, Analyzing_Institute) %>%
#   mutate(sum_count = sum(count,na.rm=TRUE)) %>%
#   mutate(sum_reads = sum(total_reads,na.rm=TRUE)) %>%
#   ungroup() %>%
#   distinct(phylum, Analyzing_Institute, .keep_all = TRUE) %>%
#   arrange(phylum, Analyzing_Institute)

