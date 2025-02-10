#!/usr/bin/env bash
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Re-blast ASVs to get latest taxonomy. Edit this file depending on the project.
# Will run MEGAN too
# 
# created May 2023
# K. Pitz MBARI
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

################################################################################
# BLAST CLUSTERS with mbon-master spawning VMs
################################################################################
# This script uses GNU parallel to spawn several blast processes on different servers.
# After the processing has been completed, the worker will copy the file to the
# /opt/blast/Worker_Output directory on the mbon-master system, and the contents of this folder will
# be moved to /Blast_VM_output

# Parameters
BLAST_INPUT="GLOMICON_seq_merged_unique.fasta"
SCRIPT_DIR="/home/mbonteam/dev"
BLAST_DB='/MBON/blastdb/nt/nt'
STANDARD_PREFIX="GLOMICON"
DB_NAME="nt"
PERCENT_IDENTITY="80"
WORD_SIZE="20"
EVALUE="1e-5"
MAXIMUM_MATCHES="100"
bitscore_sp=200
per_ID_sp=97
bitscore_gn=150
per_ID_gn=95
MINIMUM_SUPPORT="1"
TOP_PERCENT="2"
MINIMUM_SCORE="100"
LCA_PERCENT="80"
MAX_EXPECTED="1e-25"
NOTIFY_EMAIL="YES"
EMAIL_ADDRESS="kpitz@mbari.org"
 
# Define a variable called START_TIME
START_TIME=$(date +%Y%m%d_%H%M)
# make an analysis directory with starting time timestamp
# ANALYSIS_DIR="${ANALYSIS_DIRECTORY}"/Analysis_"${START_TIME}"
ANALYSIS_DIR=Analysis_"${START_TIME}"
mkdir "${ANALYSIS_DIR}"



current_time=$(date +"%T")
echo "Blast started : $current_time"
echo "Number of ASVs:"
(grep '>' "$BLAST_INPUT" | wc -l)


blast_output="${STANDARD_PREFIX}_${DB_NAME}.xml"

echo ":::::::::::::::::::::: BLAST ::::::::::::::::::::::::::::::::"
echo $(date +%Y%m%d_%H%M) " BLASTing..."
echo " Database: " "${BLAST_DB}"
echo " Percent identity: " "${PERCENT_IDENTITY}"
echo " Word size: " "${WORD_SIZE}"
echo " E value: " "${EVALUE}"
echo " Maximum target matches: " "${MAXIMUM_MATCHES}"
echo " Output format: 5"
echo " Blast input: " "${BLAST_INPUT}"
echo " Blast output: " "${blast_output}"

# Temporary destination for output files
TempDestination="/opt/blast/Worker_Output"
#remove any files that may be present in this folder
rm -f ${TempDestination}/*
# Final destination for these output files after all the blasting is done
FinalDestination="Blast_VM_output"


# The number of cores for blast to use on each system
n_cores=8
# Number of virtual machines
NumWorkers=6

# Calculate the size of each of the file chunks
#cd $FileDir
# Get the number of records in the input file
NumRecords=`grep "^>" $BLAST_INPUT | wc -l`
echo "Number of Records:"
echo "$NumRecords"
# Get the size of each chunk depending on the number of workers
SizePerWorker=$(($NumRecords / $NumWorkers))
echo "Worker File size:"
echo "$SizePerWorker"

# Add 1 just to avoid rounding errors
SizePerWorker=$(($SizePerWorker + 1 ))
echo "Worker final file size:"
echo "$SizePerWorker"


cat $BLAST_INPUT | \
parallel -S mbon-worker1,mbon-worker2,mbon-worker3,mbon-worker4,mbon-worker5,mbon-worker6 \
--N $SizePerWorker \
--recstart '>' \
--pipe blastn-wrapper \
    -db "$BLAST_DB" \
    -num_threads "$n_cores" \
    -perc_identity "${PERCENT_IDENTITY}" \
    -word_size "${WORD_SIZE}" \
    -evalue "${EVALUE}" \
    -max_target_seqs "${MAXIMUM_MATCHES}" \
    -outfmt 5 \
    -out "${blast_output}" \
    -gapopen 5 \
    -gapextend 2 \
    -reward 2 \
    -penalty -3 \
    -query -



#output should be saved in analysis directory with different appended names for worker VMs
#Will have multiple xml files as output.

# Make sure the output directory exists
mkdir $FinalDestination > /dev/null 2>&1
# Move the output files from the temporary destination to the final destination
mv ${TempDestination}/* $FinalDestination

echo $(date +%Y%m%d_%H%M) " Blast finished."

# check for blast output directory and make sure it isn't empty
[ "$(ls -A "${FinalDestination}")" ] && echo "Blast Not Empty" || echo "Blast Empty"

#list blast files that were created
(ls -ltr ${FinalDestination})

#list number of unique query IDs that were blasted, make sure it matches number of ASVs
(grep '<Iteration_query-def>' ${FinalDestination}/*.xml | wc -l)


FINISH_TIME=$(date +%Y%m%d_%H%M)

if [ "$NOTIFY_EMAIL" = "YES" ]; then
	echo 'Blast finished! Started at' $START_TIME 'and finished at' $FINISH_TIME 'Data here:' "${ANALYSIS_DIR}" | mail -s "banzai is finished" "${EMAIL_ADDRESS}"
else
	echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME
fi
echo 'Data are in ' "${ANALYSIS_DIR}"
