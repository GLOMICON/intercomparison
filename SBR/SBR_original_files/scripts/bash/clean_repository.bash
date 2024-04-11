#!/bin/bash

source ${CONFIG}

mkdir -p tmp_${PROJECT}
mv slurm* tmp_${PROJECT}
mv tmp.* tmp_${PROJECT}
mv chunck_files.txt tmp_${PROJECT}
mv asvs.fasta tmp_${PROJECT}
