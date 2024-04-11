#!/bin/bash

source $CONFIG

# One job per fastq file
sbatch --wait \
	--array 1-$(cat ${MANIFEST_DADA2_FILTERING} | wc -l)%30 \
	--export=CONFIG=${CONFIG} \
	scripts/analyses/03_dada_filtering_b.bash

sleep 1

# Gather logs
echo "sample \
in_reads \
out_reads" | tr ' ' '\t' > log/${PROJECT}_dada_filtering.log

cat log/*_${PROJECT}_round*_dada_filtering.log | sort -k2nr >> log/${PROJECT}_dada_filtering.log

 
# files for denoising 
awk '
	{gsub(/^[^_]+_|_trimmed.fastq.gz/,"",$1)}
	NR>1 && !seen[$1]++ {print $1}
	' log/${PROJECT}_dada_filtering.log > manifests/${PROJECT}_dada_denoising.txt
