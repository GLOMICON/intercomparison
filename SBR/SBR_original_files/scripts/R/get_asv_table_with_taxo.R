library(data.table)

args <- commandArgs(TRUE)

input_asv <- args[1]
input_taxo <- args[2]
output <- sub("table.tsv.gz$","table.with.taxo.idtaxa.tsv.gz",input_asv)

asv_table <- fread(input_asv)
taxo <- fread(input_taxo,
              col.names=c("amplicon","taxonomy","confidence"),
              header=FALSE)

asv_table <- merge(taxo,asv_table,by="amplicon",all.y = TRUE)[order(total,decreasing = TRUE)]

fwrite(asv_table,
       file = output,
       quote = FALSE,
       sep="\t",
       na = NA)
