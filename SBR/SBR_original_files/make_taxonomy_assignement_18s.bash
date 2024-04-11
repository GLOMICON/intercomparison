#!/bin/bash
#SBATCH --partition long
#SBATCH -o slurm-%N-%j.out
#SBATCH -e slurm-%N-%j.err

############################################################################
# script name   : make_taxonomy_assignement_18s.bash
# description   : This script annotate an ASV table using idtaxa and vsearch
# usage         : sbatch make_taxonomy_assignement_18s.bash but better to run 
#                 the lines one by one and check manualy the outputs
# authors       : Charlotte Berthelier
# contact       : charlotte.berthelier@sb-roscoff.fr
############################################################################

# Clean repository
rm -f slurm*
rm -f tmp_18s/*
rm -f tmp/*
rm -f outputs/taxo_assignment/*

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
	tmp_18s \
	manifests


# Config file
CONFIG="config/18s_ASV_table_taxonomic_assignment.cfg"

############################################
# Taxonomic assignment
############################################

# taxonomic assignment using idtaxa
job_taxo_idtaxa=$(sbatch \
        --export=CONFIG=${CONFIG} \
	--wait \
        scripts/analyses/06_taxo_assign_idtaxa_a.bash | awk '{print $4}')

# taxonomic assignment using vsearch
job_taxo_vsearch=$(sbatch --dependency=afterany:$job_taxo_idtaxa \
        --export=CONFIG=${CONFIG} \
        --wait \
        scripts/analyses/06bis_taxo_assign_vsearch_a.bash | awk '{print $4}')

job_taxo_vsearch_merging=$(sbatch --dependency=afterany:$job_taxo_vsearch \
       --export=CONFIG=${CONFIG} \
       --wait \
       scripts/analyses/06bis_taxo_assign_vsearch_c.bash | awk '{print $4}')

# Clean working space
job_cleaning=$(sbatch --dependency=afterany:$job_taxo_vsearch_merging \
	--export=CONFIG=${CONFIG} \
	--wait \
	scripts/bash/clean_repository.bash)
