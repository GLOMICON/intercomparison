#!/bin/bash

source $CONFIG

SPLE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $MANIFEST_DADA2_FILTERING)

Rscript scripts/R/dada_filtering.R $SPLE $TRUNC_R1 $TRUNC_R2 $TRUNC_Q $MAX_EE
