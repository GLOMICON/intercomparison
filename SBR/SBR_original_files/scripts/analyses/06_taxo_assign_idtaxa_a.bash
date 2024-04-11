#!/bin/bash
#
#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 30GB

source $CONFIG

# tmp dir for chunks
mkdir -p tmp/chunks
rm -f tmp/chunks/*
OUTPUTDIR="./tmp/chunks"

# Split the big fasta file into several chunks
# for parallel jobs
zcat -f $DADA2_FASTA | split -l 500 --filter='gzip > $FILE.fas.gz' - tmp/chunks/file

# List the chunk files
find ./tmp/chunks -name "*.fas.gz" -type f ! -empty > tmp/chunk_files.txt

# IDTAXA

## One job per chunk is launched
sbatch --wait \
	--array 1-$(cat tmp/chunk_files.txt | wc -l)%50 \
	--export=CONFIG=${CONFIG} \
	--cpus-per-task 4 \
	scripts/analyses/06_taxo_assign_idtaxa_b.bash

sleep 1

# Gather files
find ${OUTPUTDIR} \
        -name "*.idtaxa.out.gz" \
        -type f \
        -exec cat {} + > ${IDTAXA_OUTPUT}

# rm -rf tmp/*

# Add taxo to asv_table
if [ "$EXTRA_ANNOT" = "yes" ]
then
	Rscript scripts/R/get_asv_table_with_taxo.R $DADA2_ASV_TABLE $IDTAXA_OUTPUT
fi
