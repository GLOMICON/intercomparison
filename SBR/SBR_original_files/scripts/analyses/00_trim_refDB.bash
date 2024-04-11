#!/bin/bash
#SBATCH --cpus-per-task 1                               # number of CPUs required per task
#SBATCH --mem 4GB                                       # memory per processor
#SBATCH -p fast                                         # partition
#SBATCH -o slurm.%N.%j.out                              # STDOUT file with the Node name and the Job ID (%j= jobid)
#SBATCH -e slurm.%N.%j.err                              # STDERR file with the Node name and the Job ID


module load cutadapt/1.8.3

REFDB_DIR=archives/refdb/
INPUT=$REFDB_DIR"cyanoMarks_u1.16_RefDB_16S_pr2_taxo.fna"
PRIMERS_SET="16SV4V5"
#PRIMERS_SET="18SV4"
#PRIMERS_SET="petB"

# Define variables and output files
OUTPUT=$(echo $INPUT | sed 's/.fna/_trimmed_/')$PRIMERS_SET".fna"
LOG="log/"$(basename $INPUT | sed 's/.fna/_trimming.log/')

PRIMER_F="GTGYCAGCMGCCGCGGTAA"
PRIMER_R="CCGYCAATTYMTTTRAGTTT"
ANTI_PRIMER_R="AAACTYAAAKRAATTGRCGG" # reverse complement
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))
CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH}"

# Trim forward & reverse primers, format
cat "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
     ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
     ${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" | \
     sed '/^>/ s/;/|/g' > "${OUTPUT}"
