#!/bin/bash
#
#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB

source ${CONFIG}

# One job per fastq file
sbatch --wait --array 1-$(cat ${MANIFEST_TRIMMING} | wc -l)%30 \
	--export=CONFIG=${CONFIG} \
        scripts/analyses/01_primer_trimming_b.bash

sleep 1

# Gather logs
echo "sample status in_reads in_bp too_short too_long \
too_many_n out_reads w/adapters qualtrim_bp \
out_bp w/adapters2 qualtrim2_bp out2_bp" | \
        tr ' ' '\t' > log/${PROJECT}_primer_trimming.log

cat log/*_${PROJECT}_primer_trimming.log >> log/${PROJECT}_primer_trimming.log

rm -f log/*_${PROJECT}_primer_trimming.log

# manifest for filtering, keep only files with more than 0 reads
# if reads were not mix-orientated, only round1 should be selected

awk '$8>0 && NR>1 {print $1}' log/${PROJECT}_primer_trimming.log > manifests/${PROJECT}_dada_filtering.txt
