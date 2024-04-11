#!/bin/bash
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB

source ${1}
 # Launch using sbatch scripts/analyses/00_make_refDB.rds.bash config/make_refDB_rds_petB.cfg

Rscript scripts/R/make_refDB_rds.R $PROJECT $REFDBFNA $REFDBRDS

