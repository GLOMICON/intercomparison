library(data.table)
library(magrittr)
library(seqinr)
library(plyr)

# Rscript --verbose scripts/R/taxo_assign_vsearch.R 
args <- commandArgs(TRUE)
cyanomarks_file <- args[2]
output <- args[3]

taxuniq <- function(Y){
  Y <- unique(Y)
  if(length(Y)>1){
    sub("(^.\\:).+$","\\1*",Y[1])
  }else{
    Y
  }
}

## cyanomarks
references <- read.fasta(cyanomarks_file)
references <- getAnnot(references)

tmp <- lapply(1:length(references), function(i) {
  gsub("^>|[A-Z]+\\|tax=k\\:","",references[[i]]) %>%
    strsplit(" ") %>% unlist
})

references <- matrix(unlist(tmp), ncol=2, byrow=TRUE) %>%
  data.table

setnames(references,c("id","lineage"))

###### 
annot <- fread(args[1],header=F)
setnames(annot,c("amplicon","identity","id"))
corresp <- fread(args[5])
annot <- merge(annot,corresp,by.x="amplicon",by.y="asv")
annot <- merge(references,annot,by="id")

annot <- annot[,.(ids=paste(id,collapse = ","),
        taxo=paste(lapply(tstrsplit(lineage,","),taxuniq),collapse = ",")),
     by=list(amplicon,identity,sequence)]

fwrite(annot,paste0(output,"_",args[4]),sep="\t",row.names=F,col.names = F,quote = F,na="NA")
