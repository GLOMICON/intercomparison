library(dada2)
library(data.table)

# train idtaxa on refDB.fasta to generate refdb.rds

args <- commandArgs(TRUE)
project <- args[1]
refdbfna <- args[2]
refdbrds <- args[3]

# import refs
refs_fasta <- Biostrings::readDNAStringSet(refdbfna)

# extract taxonomy
taxo <- sub("^[^ ]+\\s","",names(refs_fasta))
taxo <- sub("^","Root;",taxo)


# train idtaxa
mzg_decipher_trained <- DECIPHER::LearnTaxa(refs_fasta, taxo)

# export the idtaxa trained refdb
saveRDS(mzg_decipher_trained, file = refdbrds)

