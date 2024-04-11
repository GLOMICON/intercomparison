#!/bin/bash

TASK_ID=$1
source $2

read SAMPLE RUNID FORWARD REVERSE <<< $(awk "NR==$TASK_ID" ${MANIFEST_TRIMMING})

ROUND_F1="outputs/reads/fwd/${SAMPLE}_${RUNID}_${PROJECT}_round1_trimmed.fastq.gz"
ROUND_F2="outputs/reads/fwd/${SAMPLE}_${RUNID}_${PROJECT}_round2_trimmed.fastq.gz"
ROUND_R1="outputs/reads/rev/${SAMPLE}_${RUNID}_${PROJECT}_round1_trimmed.fastq.gz"
ROUND_R2="outputs/reads/rev/${SAMPLE}_${RUNID}_${PROJECT}_round2_trimmed.fastq.gz"
LOG="log/${SAMPLE}_${RUNID}_${PROJECT}_primer_trimming.log"
TMP_LOG=$(mktemp --tmpdir="tmp/")

CUTADAPT="cutadapt \
	-g ${PRIMER_F} \
	-G ${PRIMER_R} \
	--report=minimal \
	--discard-untrimmed \
	--minimum-length ${MIN_LENGTH} \
	--no-indels"

$CUTADAPT \
	-o ${ROUND_F1} \
	-p ${ROUND_R1} \
	${FORWARD} \
	${REVERSE} 1> ${TMP_LOG}

awk -v a="${SAMPLE}_${RUNID}_${PROJECT}_round1" \
	'BEGIN {OFS="\t"}; NR==2{print a,$0}' \
	${TMP_LOG} > ${LOG}

$CUTADAPT \
	-o ${ROUND_F2} \
	-p ${ROUND_R2} \
	${REVERSE} \
	${FORWARD}  1> ${TMP_LOG}

awk -v a="${SAMPLE}_${RUNID}_${PROJECT}_round2" \
	'BEGIN {OFS="\t"}; NR==2{print a,$0}' \
	${TMP_LOG}  >> ${LOG}

rm -f ${TMP_LOG}
