#!/bin/bash
#
#SBATCH --mem 10GB

source ${CONFIG}

Rscript scripts/R/quality_plots.R $PROJECT $MANIFEST_QUALITY
