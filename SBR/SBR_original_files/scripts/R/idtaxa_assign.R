library(data.table)

args <- commandArgs(TRUE)

input <- args[1]
threads <- args[2] |> as.numeric()
trainingSet.path <- args[3]
output <- args[4]

trainingSet <- readRDS(trainingSet.path)

otus <- Biostrings::readDNAStringSet(input)

otus_assign <- DECIPHER::IdTaxa(otus,
		      trainingSet,
		      strand = "top",
		      threshold=50,
		      processors = threads)

taxo <- vapply(otus_assign,
	       function(X) paste(X$taxon,collapse="|"),
	       character(1))

confidence <- vapply(otus_assign,
		     function(X) paste(round(X$confidence,digits=1),collapse=";"),
		     character(1))

data.table(names(taxo),taxo,confidence) |> fwrite(file=output,
						  sep="\t",
						  quote=FALSE,
						  col.names=FALSE)

