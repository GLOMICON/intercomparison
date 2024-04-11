#!/bin/bash
#SBATCH --mem 10GB

source ${CONFIG}

Rscript scripts/R/dada_table.R $PROJECT $VERSION $MINIMUM_NUMBER_ASV $MINIMUM_STATION_NB

sleep 1

gzip "outputs/asv_table/"$PROJECT"_dada2_"$VERSION".fasta"
