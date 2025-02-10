#Begin Post_Blast_processing custom script
#Katie Pitz 021025
#Goal to run megan on blast xml results and get a taxonomy table

#manually enter fields
ANALYSIS_DIR="/home/mbonteam/MBARI/kpitz/custom_blast/GLOMICON"
MINIMUM_SCORE="100"
MAX_EXPECTED="1e-25"
TOP_PERCENT="2"
MINIMUM_SUPPORT="1"
LCA_PERCENT="80"

################################################################################
# RUN MEGAN
################################################################################

echo /usr/local/megan/tools/blast2rma --in "${ANALYSIS_DIR}/Blast_VM_output/*.xml" --format BlastXML --blastMode BlastN --out "${ANALYSIS_DIR}/Blast_VM_output"/ --minScore "${MINIMUM_SCORE}" --maxExpected "${MAX_EXPECTED}" --topPercent "${TOP_PERCENT}" --minSupport "${MINIMUM_SUPPORT}" --lcaAlgorithm naive --lcaCoveragePercent "${LCA_PERCENT}"
/usr/local/megan/tools/blast2rma --in "${ANALYSIS_DIR}"/Blast_VM_output/*.xml --format BlastXML --blastMode BlastN --out "${ANALYSIS_DIR}/Blast_VM_output"/ --minScore "${MINIMUM_SCORE}" --maxExpected "${MAX_EXPECTED}" --topPercent "${TOP_PERCENT}" --minSupport "${MINIMUM_SUPPORT}" --lcaAlgorithm naive --lcaCoveragePercent "${LCA_PERCENT}"

#Now there should be a list of .rma6 files in "${ANALYSIS_DIR}/
echo 'Now for each Megan6 file, extract the taxonomy:'
for file in "${ANALYSIS_DIR}"/*.rma6
do
  echo "$file"
	echo /usr/local/megan/tools/rma2info --in "${file}"  --read2class Taxonomy --paths --majorRanksOnly > "${file}"_tpath.txt
	/usr/local/megan/tools/rma2info --in "${file}"  --read2class Taxonomy --paths --majorRanksOnly > "${file}"_tpath.txt
done


#join files together; get all taxonomy results
cat "${ANALYSIS_DIR}"/*_tpath.txt > "${ANALYSIS_DIR}"/tpath.txt



echo 'created tpath.txt file from MEGAN results'
