#!/bin/bash
#SBATCH --partition long

source ${CONFIG}

#NRUNS=$(cat ${MANIFEST_DADA2_DENOISING} | wc -l)
NRUNS=1

if [ "$NRUNS" -gt 1 ]
then
	sbatch --array 1-$(cat ${MANIFEST_DADA2_DENOISING} | wc -l)%30 \
		--wait \
		--cpus-per-task ${THREADS} \
		--export=CONFIG=${CONFIG} \
		--mem ${MEM} \
		scripts/analyses/04_dada_denoising_b.bash
else
	sbatch --array 1 \
		--wait \
		--cpus-per-task ${THREADS} \
		--export=CONFIG=${CONFIG} \
		--mem ${MEM} \
		scripts/analyses/04_dada_denoising_b.bash
fi

sleep 1

echo "\
sample denoisedF.read denoisedR.read merged.read \
denoisedF.seq denoisedR.seq merged.seq nochim.read nochim.seq\
" | tr ' ' '\t' > log/${PROJECT}_dada_denoising.log

cat log/*_${PROJECT}_round*_dada_denoising.log | sort -k1 >> log/${PROJECT}_dada_denoising.log
rm -f log/*_${PROJECT}_round*_dada_denoising.log
