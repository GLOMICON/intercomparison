#!/bin/bash
#
#SBATCH --partition fast

source $CONFIG

RUN=$(awk "NR==$SLURM_ARRAY_TASK_ID" $MANIFEST_DADA2_DENOISING)

echo "Rscript scripts/R/dada_denoising.R $RUN $THREADS $CONCAT"

Rscript scripts/R/dada_denoising.R $RUN $THREADS $CONCAT
