BEGIN{
	FS="_"
}

/fastq\.gz/ && /^cj/{
	r1file="/shared/projects/glomicon/archives/R1R2/"$1"_"$2"_"$3"_R1.fastq.gz"
	r2file="/shared/projects/glomicon/archives/R1R2/"$1"_"$2"_"$3"_R2.fastq.gz"
	sample=$1"_"$2
	if (!a[sample]++) print sample, "Glomicon_18S", r1file, r2file
	}
	
# gsub(/-/,"_",sample)

#BEGIN{
#    FS="_"
#}
#
#/fastq\.gz/ && /^TON/{
#    r1file="archives/raw/18S/"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_18S_R1.fastq.gz"
#    r2file="archives/raw/18S/"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_18S_R2.fastq.gz"
#    sample=$2"-"$3"-"$6
#    if (!a[sample]++) print sample, "TONGA_18S", r1file, r2file
#    }
