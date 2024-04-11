library(data.table)
library(dada2)

########## arguments
args <- commandArgs(TRUE)
run_id <- args[1]
n_threads <- args[2] |> as.numeric()
just_cat <- args[3]

#####

getN <- function(x) sum(getUniques(x))
getA <- function(x) length(getUniques(x))

file_names <- dir("outputs/reads/fwd/dada_filtered/") |>
  grep(pattern = run_id, value = TRUE, fixed = TRUE)

sample_names <- sub("_trimmed.fastq.gz", "", file_names)

# File parsing
filtpathF <- "outputs/reads/fwd/dada_filtered"
filtpathR <- "outputs/reads/rev/dada_filtered"
filtFs <- paste(filtpathF, file_names, sep = "/")
filtRs <- paste(filtpathR, file_names, sep = "/")
names(filtFs) <- sample_names
names(filtRs) <- sample_names

set.seed(100)

# Learn forward error rates

errF <- learnErrors(filtFs, nbases = 1e8, multithread = n_threads)
pdf(paste0("outputs/error_plot/errF_", run_id, ".pdf"))
print(plotErrors(errF, nominalQ = TRUE))
dev.off()

# Learn reverse error rates

errR <- learnErrors(filtRs, nbases = 1e8, multithread = n_threads)
pdf(paste0("outputs/error_plot/errR_", run_id, ".pdf"))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

derepF <- derepFastq(filtFs)
ddF <- dada(derepF, err = errF, multithread = n_threads, pool = "pseudo")

derepR <- derepFastq(filtRs)
ddR <- dada(derepR, err=errR, multithread = n_threads, pool = "pseudo")

if (just_cat == "yes") {
  merger <- mergePairs(ddF,
                     derepF,
                     ddR,
                     derepR,
                     trimOverhang = TRUE,
                     returnRejects = TRUE)

  concat <- mergePairs(ddF,
                       derepF,
                       ddR,
                       derepR,
                       trimOverhang = TRUE,
                       justConcatenate = TRUE)

  for (x in names(merger)) {
    tmp <- concat[[x]][!merger[[x]]$accept & merger[[x]]$nmatch < 3, ]
    merger[[x]][!merger[[x]]$accept & merger[[x]]$nmatch < 3, ] <- tmp
    merger[[x]] <- merger[[x]][merger[[x]]$accept,]
  }
} else {
    merger <- mergePairs(ddF,
                     derepF,
                     ddR,
                     derepR,
                     trimOverhang = TRUE,
                     returnRejects = FALSE)
}

track <- data.table(
  sample = names(ddF),
  denoisedF.read = sapply(ddF, getN),
  denoisedR.read = sapply(ddR, getN),
  merged.read = sapply(merger, getN),
  denoisedF.seq = sapply(ddF, getA),
  denoisedR.seq = sapply(ddR, getA),
  merged.seq = sapply(merger, getA)
)

merger <- removeBimeraDenovo(merger, multithread = n_threads)
seqtab <- makeSequenceTable(merger)

saveRDS(seqtab, paste0("outputs/asv_table/asv_table_",run_id,".rds"))

seqtab <- t(seqtab) |> 
  data.table(keep.rownames=TRUE)
setnames(seqtab,"rn","sequence")

track[, nochim.read := sapply(sample, function(X){
  sum(seqtab[, get(X)])
})]

track[, nochim.seq := sapply(sample, function(X){
  sum(seqtab[, get(X)] != 0)
})]

fwrite(track,
       file = paste0("log/", run_id, "_dada_denoising.log"),
       col.names = FALSE,
       sep = "\t",
       quote = FALSE)
