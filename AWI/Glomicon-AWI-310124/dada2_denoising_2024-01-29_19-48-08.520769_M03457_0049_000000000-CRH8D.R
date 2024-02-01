#set software pathes
#overwrite personal library pathes: definition of system library pathes and common project library path
.libPaths(c('/isibhv/projects/p_bioinf2/Tools/common_libs/R/4.3','/usr/local/lib/R/site-library','/usr/lib/R/site-library','/usr/lib/R/library'))


# load needed packages
library(dada2)
packageVersion("dada2")
library(stringr)
packageVersion("stringr")
library(ggplot2)
packageVersion("ggplot2")
library(ShortRead)
packageVersion("ShortRead")



#get timestamp to distinguish between analysis runs
time.stamp <- gsub(" ","_",gsub(":","-",Sys.time()))

#is it a test_run? if TRUE only subset of files will be analysed (TRUE or FALSE)
test_run=FALSE

#PRIMER
#set primer sequences
FWD_PRIMER="GCGGTAATTCCAGCTCCAA"
REV_PRIMER="ACTTTCGTTCTTGAT"
#create reverse complements of primer sequences
FWD_PRIMER.RC <- dada2:::rc(FWD_PRIMER)
REV_PRIMER.RC <- dada2:::rc(REV_PRIMER)

#define project directory
project_dir <- file.path("/isibhv/projects/p_bioinf2/GLOMICON-INTERCOMP_R-4.3.2")

#set sequencing run to denoise (location of sequencing run in working directory) and raw sequence directory
seq_run <- "M03457_0049_000000000-CRH8D"
raw_dir <- file.path("/isibhv/projects/p_bioinf2/GLOMICON-INTERCOMP_R-4.3.2/raw_data",seq_run)




#PATHES AND DIRECTORIES FOR INTERMEDIATE RESULTS
#construct needed pathes
#work_dir <- file.path(project_dir,time.stamp)
temp_dir <- file.path(project_dir,paste(time.stamp,seq_run,"temp",sep="_"))
preFilt_dir <- file.path(temp_dir,"preFilt")
primerCut5_dir <- file.path(temp_dir,"primerCut5")
primerCut3_dir <- file.path(temp_dir,"primerCut3")
qualFiltTrim_dir <- file.path(temp_dir,"qualFiltTrim")
#create directories
#if(!dir.exists(work_dir)) dir.create(work_dir)
if(!dir.exists(temp_dir)) dir.create(temp_dir)
if(!dir.exists(preFilt_dir)) dir.create(preFilt_dir)
if(!dir.exists(primerCut5_dir)) dir.create(primerCut5_dir)
if(!dir.exists(primerCut3_dir)) dir.create(primerCut3_dir)
if(!dir.exists(qualFiltTrim_dir)) dir.create(qualFiltTrim_dir)


#construct needed file lists
#basenames of all files
fnFs.files <- basename(sort(list.files(raw_dir, pattern ="_L001_R1_001.fastq.gz", full.names = TRUE)))
fnRs.files <- basename(sort(list.files(raw_dir, pattern ="_L001_R2_001.fastq.gz", full.names = TRUE)))
#basenames of a file subset
fnFs.files.sub <- sample(fnFs.files,4)
fnRs.files.sub <- gsub("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz",fnFs.files.sub, fixed=TRUE)
#pathes of raw files to analyse
fnFs.raw <- file.path(raw_dir,fnFs.files)
fnRs.raw <- file.path(raw_dir,fnRs.files)
#use sub set of file pathes if test_run is TRUE
if (test_run==TRUE)
  {
  fnFs.raw <- file.path(raw_dir,fnFs.files.sub)
  fnRs.raw <- file.path(raw_dir,fnRs.files.sub)
  }
fnFs.preFilt <- file.path(preFilt_dir,basename(fnFs.raw))
fnRs.preFilt <- file.path(preFilt_dir,basename(fnRs.raw))
fnFs.primerCut5 <- file.path(primerCut5_dir,basename(fnFs.raw))
fnRs.primerCut5 <- file.path(primerCut5_dir,basename(fnRs.raw))
fnFs.primerCut3 <- file.path(primerCut3_dir,basename(fnFs.raw))
fnRs.primerCut3 <- file.path(primerCut3_dir,basename(fnRs.raw))
fnFs.qualFiltTrim <- file.path(qualFiltTrim_dir,basename(fnFs.raw))
fnRs.qualFiltTrim <- file.path(qualFiltTrim_dir,basename(fnRs.raw))


#create sample names
sample.names <- unname(sapply(fnFs.files,function(filenames) 
  paste(strsplit(str_remove(filenames,"_L001_R1_001.fastq.gz"),"_")[[1]][c(2,3,6)], collapse="_")))
sample.names
#example:
#GLOMICON-INTERCOMP_Framstrait_13_AWI_M03457_0008_000000000-C27WW_RV-AWI13_S56_L001_R1_001.fastq.gz
#Project_location_sample-number_research-institute_sequencer-ID_run-Index_flowcell-ID_sample-ID_SampleSheet-Index_lanw-number_read-index_suffix
#take substring 2 and 3 (location and sample number)

#create file to sample name mapping
file.list <- cbind(sample.names,fnFs.files, fnRs.files)
file.list

# create sample names of subsampled file set
sample.names.sub <- unname(sapply(fnFs.files.sub,function(filenames) 
  paste(strsplit(str_remove(filenames,"_L001_R1_001.fastq.gz"),"_")[[1]][c(2,3,6)], collapse="_")))
sample.names.sub
#create file to sample name mapping of subsampled file list
#fnFs.files.sub <- basename(fnFs.sub)
#fnRs.files.sub <- basename(fnRs.sub)
file.list.sub <- cbind(sample.names.sub,fnFs.files.sub, fnRs.files.sub)
file.list.sub

# if test run change file list and sample names  by their subsampled versions
if (test_run==TRUE)
  {
  file.list <- file.list.sub
  sample.names <- sample.names.sub
  }

#QUALITY PLOTS OF RAW FILES
#plot Phred quality score profiles of subset of raw files; data are aggregated; always subset to reduce computational costs
#plot R1 data
plotQualityProfile(file.path(raw_dir,fnFs.files.sub),aggregate=TRUE)
#plot R2 data
plotQualityProfile(file.path(raw_dir,fnRs.files.sub),aggregate=TRUE)


#PREFILTERING
prefilterOut <- filterAndTrim(fnFs.raw,fnFs.preFilt,fnRs.raw,fnRs.preFilt,truncQ=2,minLen=150,maxN=0,multithread = TRUE)
prefilterOut

#PREPARE primer checks; create all possible primer versions (fwd, fwd_complement, rev, rev_complement)
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  #require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  
}

FWD_PRIMER.orients <- allOrients(FWD_PRIMER)
REV_PRIMER.orients <- allOrients(REV_PRIMER)

FWD_PRIMER.orients 
REV_PRIMER.orients 


#check primer situation in sample subset of sequence files

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

for(i in 1:nrow(file.list.sub)){
  colnames(file.list.sub) <- NULL
  print(file.list.sub[i,1])
  print(rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = file.path(preFilt_dir,file.list.sub[i,2])), 
        FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = file.path(preFilt_dir,file.list.sub[i,3])), 
        REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = file.path(preFilt_dir,file.list.sub[i,2])), 
        REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = file.path(preFilt_dir,file.list.sub[i,3]))))
}


#REMOVE PRIMERS at 5'-end with cutadapt
#which cutadapt version?
system2("cutadapt", args = "--version")

#primer sequences
FWD_PRIMER
REV_PRIMER

# Run Cutadapt to trim both primer at 5' end of R1 and R2 respectively
for(i in seq_along(fnFs.preFilt)) {
  system2("cutadapt", args = c("-g", paste("\"",FWD_PRIMER,";min_overlap=",str_length(FWD_PRIMER),";max_error_rate=0.2","\"",sep=""), 
                             "-G", paste("\"",REV_PRIMER,";min_overlap=",str_length(REV_PRIMER),";max_error_rate=0.2","\"",sep=""),
                             "--discard-untrimmed", "--cores=0", "--minimum-length 150","-o", fnFs.primerCut5[i], "-p", fnRs.primerCut5[i],
                             fnFs.preFilt[i], fnRs.preFilt[i]))
  
}


#REMOVE primer reverse complements at 3'-end (affected sequences might be artifacts anyway because they are too short, depending on target)

#reverse complements of primer sequences
FWD_PRIMER.RC
REV_PRIMER.RC

# Run Cutadapt both primer at 3' end
for(i in seq_along(fnFs.primerCut5)) {
  system2("cutadapt", args = c("-a", paste("\"",REV_PRIMER.RC,";min_overlap=",str_length(REV_PRIMER.RC),";max_error_rate=0.2","\"",sep=""), 
                             "-A", paste("\"",FWD_PRIMER.RC,";min_overlap=",str_length(FWD_PRIMER.RC),";max_error_rate=0.2","\"",sep=""),
                             "--cores=0", "--minimum-length 150", "-o", fnFs.primerCut3[i], "-p", fnRs.primerCut3[i],
                             fnFs.primerCut5[i], fnRs.primerCut5[i]))
  
}




#check primer situation in sample subset of trimmed sequences files
FWD_PRIMER.orients 
REV_PRIMER.orients 

for(i in 1:nrow(file.list.sub)){
  colnames(file.list.sub) <- NULL
  print(file.list.sub[i,1])
  print(rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = file.path(primerCut3_dir,file.list.sub[i,2])), 
              FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = file.path(primerCut3_dir,file.list.sub[i,3])), 
              REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = file.path(primerCut3_dir,file.list.sub[i,2])), 
              REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = file.path(primerCut3_dir,file.list.sub[i,3]))))
}





#plot Phred quality score profiles of subset of trimmed sequence files; data are aggregated
#plot R1 data
plotQualityProfile(file.path(primerCut3_dir,fnFs.files.sub),aggregate=TRUE)
#plot R2 data
plotQualityProfile(file.path(primerCut3_dir,fnRs.files.sub),aggregate=TRUE)


#filter and trim sequences; truncLen parameters estimated by inspecting quality profiles
#roughly when green line (phred score mean) drops below 30
#maxEE is considered to be a function of truncLen; it is easiest to set a factor/proportion between both values, e.g. 100
filterOut <- filterAndTrim(fnFs.primerCut3,fnFs.qualFiltTrim,
                           fnRs.primerCut3,fnRs.qualFiltTrim, 
                           maxN=0,maxEE=c(2.70,2.30),
                           truncLen=c(270,230),
                           verbose = TRUE, rm.phix = TRUE, 
                           compress = TRUE, multithread = TRUE)

# sequences lost by quality filtering
filterOut


#plot Phred quality score profiles of subset of quality filtered sequence files; data are aggregated
#plot R1 data
plotQualityProfile(file.path(qualFiltTrim_dir,fnFs.files.sub),aggregate=TRUE)
#plot R2 data
plotQualityProfile(file.path(qualFiltTrim_dir,fnRs.files.sub),aggregate=TRUE)



# dereplicate identical sequences to reduce computational costs
fnFs.deRep <- derepFastq(fnFs.qualFiltTrim)
fnRs.deRep <- derepFastq(fnRs.qualFiltTrim)
names(fnFs.deRep) <- sample.names
names(fnRs.deRep) <- sample.names



#LEARN DADA2 ERROR RATES
errFWDs <- learnErrors(fnFs.deRep, multithread=TRUE,randomize=TRUE, nbases = 1e9)
errREVs <- learnErrors(fnRs.deRep, multithread=TRUE,randomize=TRUE, nbases = 1e9)

#plot error profiles
#black curve should fit black dots and monotonically decreasing
plotErrors(errFWDs, nominalQ=TRUE)
plotErrors(errREVs, nominalQ=TRUE)


#DADA SAMPLE INFERENCE
#pool="pseudo" if computation too costly
dadaFWDs <- dada(fnFs.deRep, err=errFWDs, multithread=TRUE, pool=TRUE)
dadaREVs <- dada(fnRs.deRep, err=errREVs, multithread=TRUE, pool=TRUE)

#inspect dada inference objects
dadaFWDs
dadaREVs

#MERGE PAIRED ENDS
mergers <- mergePairs(dadaFWDs, fnFs.deRep, dadaREVs, fnRs.deRep, minOverlap=40, maxMismatch=0, verbose=TRUE)
str(mergers)

#CONSTRUCT, INSPECT AND SAVE SEQUENCE TABLE
#create seqtab
seqtab <- makeSequenceTable(mergers)
#number of samples
dim(seqtab)[1]
#number of ASVs
dim(seqtab)[2]
#length distribution of ASVs
table(nchar(getSequences(seqtab)))
#save sequence tabs
saveRDS(seqtab, file = paste(time.stamp,"_",seq_run,"_seqtab",".rda",sep=""))
write.csv(t(seqtab), paste(time.stamp,"_",seq_run,"_seqtab",".csv",sep=""), quote=FALSE )



#TRACK READS
getN <- function(x) sum(getUniques(x))
track <- cbind(prefilterOut,filterOut, sapply(dadaFWDs, getN), sapply(dadaREVs, getN), sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("raw","preFiltered","primerFiltered", "qualFiltered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
track
saveRDS(track, file = paste(time.stamp,"_",seq_run,"_trackedReads",".rda",sep=""))
write.csv(track, paste(time.stamp,"_",seq_run,"_trackedReads",".csv",sep=""), quote=FALSE )






