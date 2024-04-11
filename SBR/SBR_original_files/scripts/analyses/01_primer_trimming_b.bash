#!/bin/bash
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 1
#SBATCH --mem 5GB

bash scripts/bash/primer_trimming.bash $SLURM_ARRAY_TASK_ID $CONFIG

