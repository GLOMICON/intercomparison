#!/bin/bash
#SBATCH --partition long
#SBATCH -o slurm-%N-%j.out
#SBATCH -e slurm-%N-%j.err

############################################################################
# script name   : make_asv_table_18s.bash
# description   : This script produces an ASV table using dada2
# usage         : sbatch make_asv_table_18s.bash but better to run 
#                 the lines one by one and check manualy the outputs
# authors       : Charlotte Berthelier
# contact       : charlotte.berthelier@sb-roscoff.fr
############################################################################

## Clear out working directory

    rm -f manifests/*
	rm -f log/*
	rm -f outputs/asv_table/*
    rm -f outputs/asv_table/*
    rm -f outputs/error_plot/*
	rm -f outputs/reads/fwd/*
	rm -f outputs/reads/rev/*
    rm -f outputs/reads/fwd/dada_filtered/*
	rm -f outputs/reads/rev/dada_filtered/*

# Outputs, log files and
# temporary files folders are created

mkdir -p \
	archives/raw \
	archives/refdb \
	outputs/reads/fwd/dada_filtered \
	outputs/reads/rev/dada_filtered \
	outputs/taxo_assignment \
	outputs/asv_table \
	outputs/quality_plot \
	outputs/error_plot \
	log \
  	tmp \
	manifests


# Config file
CONFIG="config/18s_ASV_table_generation.cfg"


############################################
# Primer trimming
############################################
# Create the manifest

ls ../archives/R1R2/ | \
	awk -f scripts/awk/18s_primer_manifest.awk > manifests/18s_primer_trimming.txt

# Remove the primers
# parameters can be directly edited
# in config/primer_trimming_18s.cfg

job_id_trimming=$(sbatch \
	--export=CONFIG=${CONFIG} \
	--wait \
	scripts/analyses/01_primer_trimming_a.bash | awk '{print $4}') 


############################################
# Quality assesment
############################################

# Generate quality plots
job_quality_plot=$(sbatch --dependency=afterany:$job_id_trimming \
	--export=CONFIG=${CONFIG} \
	scripts/analyses/02_quality_plot.bash )


############################################
# Dada2 filtering
############################################

# run filtering
job_id_filtering=$(sbatch --dependency=afterany:$job_id_trimming \
	--export=CONFIG=${CONFIG} \
	--wait \
	scripts/analyses/03_dada_filtering_a.bash | awk '{print $4}')


############################################
# Dada2 denoising
############################################

job_id_denoising=$(sbatch --dependency=afterany:$job_id_filtering \
	--export=CONFIG=${CONFIG} \
	--wait \
	scripts/analyses/04_dada_denoising_a.bash | awk '{print $4}')

############################################
# ASV table assembling
############################################

job_id_table=$(sbatch --dependency=afterany:$job_id_denoising \
	--export=CONFIG=${CONFIG} \
	--wait \
	scripts/analyses/05_dada_table.bash | awk '{print $4}')


# Remove "NNNNNNNNNN" in sequences if CONCAT="yes"
CONCAT="yes"
if [ "$CONCAT" == "yes" ]; then
	job_clean_table=$(sbatch --dependency=afterany:$job_id_table \
	--export=CONFIG=${CONFIG} \
	scripts/bash/clean_sequences.bash)
fi

job_cleaning=$(sbatch --dependency=afterany:$job_id_table \
	--export=CONFIG=${CONFIG} \
	--wait scripts/bash/clean_repository.bash )
