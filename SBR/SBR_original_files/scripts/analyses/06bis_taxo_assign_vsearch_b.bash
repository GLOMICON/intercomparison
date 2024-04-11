#!/bin/bash

#SBATCH --partition fast
#SBATCH --cpus-per-task 1
#SBATCH --mem 10GB

source ${CONFIG}

TMPINPUT=$(mktemp --tmpdir=${my_wk})
TMPCORRESP=$(mktemp --tmpdir=${my_wk})
TMPOUTPUT=$(mktemp --tmpdir=${my_wk})
TMPOUTPUT2=$(mktemp --tmpdir=${my_wk})

if [[ $SLURM_ARRAY_TASK_ID -lt $NCHUNK ]]
then
	awk 'NR=='$((($SLURM_ARRAY_TASK_ID - 1) * $CHUNK_SIZE + 1))',NR=='$(($SLURM_ARRAY_TASK_ID * $CHUNK_SIZE))'' ${TMP_FASTA} > ${TMPINPUT}
else
	awk 'NR=='$((($SLURM_ARRAY_TASK_ID - 1) * $CHUNK_SIZE + 1))',NR=='$((($SLURM_ARRAY_TASK_ID - 1) * $CHUNK_SIZE + $NLINES_LAST_CHUNK))'' ${TMP_FASTA} > ${TMPINPUT}
fi

echo 'asv sequence' > ${TMPCORRESP}
awk 'BEGIN{RS=">";OFS=" "}NR>1{print $1,$2}' ${TMPINPUT} >> ${TMPCORRESP}

vsearch --usearch_global ${TMPINPUT} \
	-db ${REFDB_VSEARCH} \
	--id ${ID_VSEARCH} \
	--maxaccepts 0 ${TOP_HITS_ONLY} \
	--userout ${TMPOUTPUT} \
	--maxhits ${MAX_HITS} \
	--iddef ${IDDEF} \
	--userfields ${USERFIELDS}

sed 's/;size=/\t/
        s/;//' ${TMPOUTPUT} > ${TMPOUTPUT2}

Rscript --verbose /shared/projects/tonga_metab/01_dada2_asv_table_generation_and_annotation/scripts/R/taxo_assign_vsearch.R ${TMPOUTPUT2} ${REFDB_VSEARCH} ${VSEARCH_OUTPUT} ${SLURM_ARRAY_TASK_ID} ${TMPCORRESP} > log/taxo_assign_vsearch.Rout

sleep 2
