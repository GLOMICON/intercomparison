library(data.table)

args <- commandArgs(TRUE)
project <- args[1]
manifest <- args[2]

# Files to work with
input <- fread(manifest)
input[,run := gsub("^[^_]+_|_round.$", "", sample)]
input <- input[out_reads>0, ][order(sample)]
input[, files := paste0(sample, "_trimmed.fastq.gz")]

pathF <- "outputs/reads/fwd"
pathR <- "outputs/reads/rev"
path_output <- "outputs/quality_plot"

## indiv plots

pdf(file.path(path_output, paste(project, "indiv_F_Qplots.pdf", sep = "_")))
for(i in file.path(pathF, input[, files])){
  print(dada2::plotQualityProfile(i))
}
dev.off()

pdf(file.path(path_output, paste(project, "indiv_R_Qplots.pdf", sep = "_")))
for(i in file.path(pathR, input[, files])){
  print(dada2::plotQualityProfile(i))
}
dev.off()

## summarized plots

summarized_plots <- function(files,prefix){
  
  fastqs_round1 <- grep("_round1_", files, value = TRUE)
  fastqs_round2 <- grep("_round2_", files, value = TRUE)
  
  if(length(fastqs_round1)>0) {
    png(file.path(path_output, paste(prefix, "F_round1_Qplots.png", sep = "_")))
    print(dada2::plotQualityProfile(file.path(pathF, fastqs_round1),
                                    aggregate = TRUE))
    dev.off()
    
    png(file.path(path_output, paste(prefix, "R_round1_Qplots.png", sep = "_")))
    print(dada2::plotQualityProfile(file.path(pathR, fastqs_round1),
                                    aggregate = TRUE))
    dev.off()
  }
  
  if(length(fastqs_round2) > 0) {
    png(file.path(path_output, paste(prefix, "F_round2_Qplots.png", sep = "_")))
    print(dada2::plotQualityProfile(file.path(pathF, fastqs_round2),
                                    aggregate = TRUE))
    dev.off()
    
    png(file.path(path_output,paste(prefix, "R_round2_Qplots.png", sep = "_")))
    print(dada2::plotQualityProfile(file.path(pathR, fastqs_round2),
                                    aggregate = TRUE))
    dev.off()
  }
}

# all together
summarized_plots(files = input[,files], prefix = project)

# per run
lapply(split(input, by = "run"), function(X){
  summarized_plots(files = X$files,
                   prefix = paste(project, unique(X$run), sep = "_"))
})
