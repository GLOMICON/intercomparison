# taxonomic assignment with PR2 and dada2
# adapted from script by Daniel Vaulot
# KP 041525

library(dada2); packageVersion("dada2")
library(dplyr)

# example from dada2 tutorial (with our sequences)
# # list of sequences to assign
# seqs <- c('CAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATTTCGGGTCGGGCCGAGCGGTCTGCCGATGGGTATGCACTGTTTGGCGCGGCCTTCTTTCCGGAGACCGCGGCTACTCTTAACTGAGCGGGCGTGGGAGACGGATCGTTTACTTTGAAAAAATCAGAGTGTTTCTAGCAGGCAGCTCGCTCTTGCATAGGTTAGCATGGGATAATTTAATAGGACTCTGGTGCTATTTTGTTGGTTTCGAACACCGGAGTAATGATTAAAAGGGGCAGTCAGGGGCACTCGTATTCCGTCGAGAGAGGTGAAATTCTCAGACCAATGGAAGACGAACCACTGC',
#           'GCACCTACCGATTGAATGGTCCGGTGAAGACTCGGGATTGTGGTCTGGCTCCTTCATTGGGGCCAGACCGTGAGAACTTGTCTGAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC',
#           'GCACCTACCGATTGAATGGTCCGGTGAGGCCTCGGGATCGTGGCGAACTTTCTTCATTGGAGGTGAGCTGTGAGAACTTGTCCAAATCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC',
#           'GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGGTTGGTTTCCTTTATTGGAATCTGACCACGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC')
# 
# #dada2 PR2 fasta file: (not trained, just fasta file)
# file <- '/Users/kpitz/Projects/PR2/pr2_version_5.1.0_SSU_dada2.fasta'
# 
# set.seed(100) # Initialize random number generator for reproducibility
# taxa <- assignTaxonomy(seqs, file, multithread=FALSE, taxLevels = c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family","Genus","Species"))
# #assignTaxonomy(..., taxLevels = c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family","Genus","Species"))
# unname(taxa)

## From DV example --------------------

dada2_assign <- function(seq_file_name,
                         ref_file_name,
                         tax_levels){
  
  # It is necessary to read the sequences to get the names because dada2 takes the sequences themselves as names.
  seq_in <- Biostrings::readDNAStringSet(seq_file_name)
  seq_names <- names(seq_in)
  
  taxa <- dada2::assignTaxonomy(seqs=seq_file_name,
                                refFasta=ref_file_name,
                                taxLevels = tax_levels,
                                minBoot = 0, 
                                outputBootstraps = TRUE,
                                multithread = TRUE,
                                verbose = TRUE)
  
  cat("\nDone with assign\n")
  
  boot <- data.frame(taxa$boot) %>% 
    dplyr::rename_all(funs(stringr::str_c(.,"_boot")))
  
  dada2_result <- data.frame(sequence_hash = seq_names) %>%
    dplyr::bind_cols(data.frame(taxa$tax)) %>% 
    dplyr::bind_cols(boot)
  
  cat("Write file\n")
  readr::write_tsv(dada2_result, filename_change_ext(seq_file_name,"dada2.taxo"), na="")
  
  return(dada2_result)
  
}

filename_change_ext <- function (file_name, new_ext){
  
  # Do it twice for files with 2 extensions sur as .fasta.gz
  file_name_out <- fs::path_ext_set(fs::path_ext_remove(file_name), new_ext)
}

## RUN ----------------------------
tax_levels <- c("domain", "supergroup", "division", "subdivision",  
                "class", "order", "family", "genus", "species") 

# 13542 18S ASV sequences
seq_file_name <- '/home/kpitz/GLOMICON/GLOMICON_seq_merged_unique.fasta'
ref_file_name <- '/home/kpitz/ref_db/PR2/pr2_version_5.1.0_SSU_dada2.fasta'

result <- dada2_assign(seq_file, ref_file, tax_levels)

