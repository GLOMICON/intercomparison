library(data.table)

### arguments
args <- commandArgs(TRUE)
sple <- args[1]

if(grepl("_round2$",sple)){
	trunc_fwd <- args[3] |> as.numeric()
	trunc_rev <- args[2] |> as.numeric()
}else{
	trunc_fwd <- args[2] |> as.numeric()
	trunc_rev <- args[3] |> as.numeric()
}

trunc_q <- args[4] |> as.numeric()

maxee <- args[5] |> as.numeric()

fastq <- paste(sple, "trimmed.fastq.gz", sep = "_")
output_log <- paste0("log/",sple,"_dada_filtering.log")

# File parsing
pathF <- "outputs/reads/fwd"
pathR <- "outputs/reads/rev"
filtpathF <- file.path(pathF, "dada_filtered")
filtpathR <- file.path(pathR, "dada_filtered")


# Filtering
dada2::filterAndTrim(
	fwd = file.path(pathF, fastq),
	filt = file.path(filtpathF, fastq),
	rev = file.path(pathR, fastq),
	filt.rev = file.path(filtpathR, fastq),
	truncLen = c(trunc_fwd, trunc_rev),
	truncQ = trunc_q,
	maxEE = maxee,
	maxN = 0,
	compress = TRUE,
	verbose = TRUE,
	multithread = FALSE
) -> out

# export
data.table(out, keep.rownames = TRUE) |>
  fwrite(output_log, sep = "\t", col.names = FALSE)
