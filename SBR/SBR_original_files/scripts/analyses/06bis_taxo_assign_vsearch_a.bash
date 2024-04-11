#!/bin/bash
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 1
#SBATCH --mem 10GB

source $CONFIG

# Temporary files
TMP_FASTA="asvs.fasta"
#ASVS_TMP=$(mktemp --tmpdir="/shared/projects/tonga_metab/01_dada2_asv_table_generation_and_annotation/")
ASVS_TMP=$(mktemp --tmpdir=${my_wk})

## generate filtered fasta file out of the dada2 filtered.table.tsv.gz
zcat "${DADA2_ASV_TABLE_vsearch}" | \
	            awk 'NR>1 {print $2}' >> "${ASVS_TMP}"

cat "${ASVS_TMP}" | \
	sort -u | \
	awk '{print ">ASV"NR"\n"$1}' > "${TMP_FASTA}"

echo -n "" > ${VSEARCH_OUTPUT}

NLINES_FASTA=$(wc -l < $TMP_FASTA)

CHUNK_SIZE=1000


if (($NLINES_FASTA % $CHUNK_SIZE))
then
	NCHUNK=$(($NLINES_FASTA / $CHUNK_SIZE))
	NLINES_LAST_CHUNK=$CHUNK_SIZE
else
	NCHUNK=$(($NLINES_FASTA / $CHUNK_SIZE + 1))
	NLINES_LAST_CHUNK=$(($NLINES_FASTA % $CHUNK_SIZE))
fi

sbatch --array 1-${NCHUNK} \
	--constraint=intel \
	--wait \
	--export=TMP_FASTA=${TMP_FASTA},CONFIG=${CONFIG},VSEARCH_OUTPUT=${VSEARCH_OUTPUT},REFDB_VSEARCH=${REFDB_VSEARCH},NCHUNK=${NCHUNK},CHUNK_SIZE=${CHUNK_SIZE},NLINES_LAST_CHUNK=${NLINES_LAST_CHUNK} \
	scripts/analyses/06bis_taxo_assign_vsearch_b.bash
