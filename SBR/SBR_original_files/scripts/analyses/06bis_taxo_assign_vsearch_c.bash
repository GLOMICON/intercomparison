#!/bin/bash

#SBATCH --partition fast
#SBATCH --cpus-per-task 1
#SBATCH --mem 10GB

source ${CONFIG}

# mv slurm* ${TMPDIR}

echo "find . -name ${VSEARCH_OUTPUT}_*  -type f -exec cat {} + > ${VSEARCH_OUTPUT}"

find outputs/taxo_assignment/ -wholename "${VSEARCH_OUTPUT}_*" -type f -exec cat {} + > "${VSEARCH_OUTPUT}"
gzip "${VSEARCH_OUTPUT}"

mv outputs/taxo_assignment/"${VSEARCH_OUTPUT}_*" tmp_${PROJECT}
# rm -rf outputs/taxo_assignment/"${VSEARCH_OUTPUT}_*"
