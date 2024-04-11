library(dada2)
library(data.table)

########## arguments
args <- commandArgs(TRUE)
project <- args[1]
version <- args[2]
MINIMUM_NUMBER_ASV <- args[3]
MINIMUM_STATION_NB <- args[4]
output_table <- paste0("outputs/asv_table/",project,"_dada2_",version,".filtered.table.tsv.gz")
output_fasta <- paste0("outputs/asv_table/",project,"_dada2_",version,".fasta")


x <- grep(pattern = paste0("asv_table_.+_",project,"_.+\\.rds$"),
	  dir("outputs/asv_table", full.names = TRUE),
	  value = TRUE)

if(length(x) > 1) {
	asv_table <- mergeSequenceTables(tables = x)
} else {
	asv_table <- readRDS(x)
}

samples <- row.names(asv_table)

asv_table <- t(asv_table) |> data.table(keep.rownames=TRUE)
setnames(asv_table,"rn","sequence")
asv_table[,amplicon:=sapply(sequence,digest::digest,algo="sha1")]

# sum samples
final_samples <- sub("_.+$","",samples) |> unique()

for(i in final_samples){
  x <- grep(paste0(i, "_"), colnames(asv_table), value= TRUE)
  asv_table[,(i):=rowSums(.SD),.SDcols=x]
  asv_table[,(x):=NULL]
}

x <- rowSums(asv_table[,.SD,.SDcols=final_samples])
asv_table[,total:=x]
x <- apply(asv_table[,.SD,.SDcols=final_samples],1,function(X) sum(X>0))
asv_table[,spread:=x]
asv_table <- asv_table[total>=MINIMUM_NUMBER_ASV & spread>=MINIMUM_STATION_NB][order(total,decreasing = T)]
asv_table <- asv_table[,.SD,.SDcols=c("amplicon","sequence","total","spread",final_samples)]

fwrite(asv_table,output_table,sep="\t",na="NA",quote=FALSE)

seqinr::write.fasta(sequences=as.list(asv_table[,sequence]),
		    names=asv_table[,amplicon],
		    file.out=output_fasta)


