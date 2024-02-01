#PREPARE ENVIRONMENT
#define project and database directories
project_dir <- file.path("/isibhv/projects/p_bioinf2/GLOMICON-INTERCOMP_R-4.3.2")
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs/")
#set software pathes
#overwrite personal library pathes: definition of system library pathes and common project library path
.libPaths(c('/isibhv/projects/p_bioinf2/Tools/common_libs/R/4.3','/usr/local/lib/R/site-library','/usr/lib/R/site-library','/usr/lib/R/library'))

#load dada2 package
library(dada2);packageVersion("dada2")

#INPUT
#define runs as input
runID.1 <- "2024-01-29_19-48-00.748157_M03457_0008_000000000-C27WW"
runID.2 <- "2024-01-29_19-48-08.520769_M03457_0049_000000000-CRH8D"

#CREATE MERGED ASV TABLE
#read seqtab rda files 
seqtab.1 <- readRDS(file.path(project_dir, paste(runID.1,"_seqtab.rda",sep="")))
seqtab.2 <- readRDS(file.path(project_dir, paste(runID.2,"_seqtab.rda",sep="")))
#Merge multiple seqtabs
seqtab.merged <- mergeSequenceTables(seqtab.1, seqtab.2)

# REMOVE CHIMERA
seqtab.merged.nochim <- removeBimeraDenovo(seqtab.merged, method="consensus", multithread=TRUE)
#save clean ASV table
write.csv(seqtab.merged.nochim, paste(basename(project_dir),"_seqtab.merged.nochim.csv",sep=""), quote=FALSE )
saveRDS(seqtab.merged.nochim, file = paste(basename(project_dir),"_seqtab.merged.nochim.rda",sep=""))


#TRACK READS
#read trackedReads rda files from different sequencing runs as dataframes
trackedReads.1 <- readRDS(file.path(project_dir,paste(runID.1,"_trackedReads.rda",sep="")))
trackedReads.2 <- readRDS(file.path(project_dir,paste(runID.2,"_trackedReads.rda",sep="")))
#merge them
trackedReads.merged <- rbind(trackedReads.1, trackedReads.2)
#add sequence numbers after chimera removal to dataframe
tracked.reads.merged.nochim <- cbind(trackedReads.merged, rowSums(seqtab.merged.nochim))
#add column headers
colnames(tracked.reads.merged.nochim) <- c("raw", "prefiltered", "primerFiltered", "qualfiltered", "denoisedF", "denoisedR", "merged", "nonchim")
#display read track table
(tracked.reads.merged.nochim)
#save read track table
saveRDS(tracked.reads.merged.nochim, file = paste(basename(project_dir),"_trackedReads",".rda",sep=""))
write.csv(tracked.reads.merged.nochim, paste(basename(project_dir),"_trackedReads",".csv",sep=""), quote=FALSE )

#TAXONOMY
#assign Silva taxonomy
taxaSilva <- assignTaxonomy(seqtab.merged.nochim, file.path(db_dir,"silva_nr_v138_train_set_Euk.fa"), taxLevels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),multithread=TRUE)
#save Silva taxonomy
write.csv(taxaSilva, paste(basename(project_dir),"_taxaSilva.csv",sep=""), quote=FALSE )
saveRDS(taxaSilva, file = paste(basename(project_dir),"_taxaSilva.rda",sep=""))
#Remove sequence rownames (for display only)
taxaSilva.print <- taxaSilva 
rownames(taxaSilva.print) <- NULL
#show first 10 lines of annotated taxonomy
head(taxaSilva.print,10)
#save annotated ASV table
write.csv(cbind(t(seqtab.merged.nochim), taxaSilva), paste(basename(project_dir),"_seqtab.merged.nochim.Silva-138.csv",sep=""), quote=FALSE )
saveRDS(cbind(t(seqtab.merged.nochim), taxaSilva), file = paste(basename(project_dir),"_seqtab.merged.nochim.Silva-138.rda",sep=""))

#assign PR2 taxonomy
taxaPR2 <- assignTaxonomy(seqtab.merged.nochim, file.path(db_dir,"pr2_version_5.0.0_SSU_dada2.fasta"), taxLevels = c("Domain", "Kingdom", "Supergroup","Division", "Class", "Order", "Family", "Genus", "Species"),multithread=TRUE)
#save PR2 taxonomy
write.csv(taxaPR2, paste(basename(project_dir),"_taxaPR2.csv",sep=""), quote=FALSE )
saveRDS(taxaPR2, file = paste(basename(project_dir),"_taxaPR2.rda",sep=""))
#Remove sequence rownames (for display only)
taxaPR2.print <- taxaPR2
rownames(taxaPR2.print) <- NULL
#show first 10 lines of annotated taxonomy
head(taxaPR2.print,10)
#save annotated ASV table
write.csv(cbind(t(seqtab.merged.nochim), taxaPR2), paste(basename(project_dir),"_seqtab.merged.nochim.PR2-500.csv",sep=""), quote=FALSE )
saveRDS(cbind(t(seqtab.merged.nochim), taxaPR2), file = paste(basename(project_dir),"_seqtab.merged.nochim.PR2-500.rda",sep=""))



