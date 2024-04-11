#!/bin/bash
#
#SBATCH --partition fast
#SBATCH --mem 20GB

TMPINPUT=$(awk "NR==$SLURM_ARRAY_TASK_ID" tmp/chunk_files.txt)

source ${CONFIG}

OUTPUT=$(echo $TMPINPUT | \
        awk -v outputdir="./tmp/chunks" -F"/" '{sub(/.fas(ta)*(\.gz)*/,".idtaxa.out.gz",$NF)
                        print outputdir"/"$NF}')

Rscript scripts/R/idtaxa_assign.R $TMPINPUT $THREADS $REFDB_IDTAXA $OUTPUT
