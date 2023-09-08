module purge

#SETTING VARIABLES
VSEARCH=$(which vsearch)
THREADS=1
#CLUSTER=0.97
LOW_ABUND_SAMPLE=2
PROJ="Aqua_Fungi_5_8S_Sweden"
TMP_FASTA1=$(mktemp)
TMP_FASTA2=$(mktemp)
TMP_FASTA7=$(mktemp)

#rename fasta headers in single files
cat *.fas >> "${TMP_FASTA7}"
perl /cluster/projects/nn9745k/scripts/bioinfo_pipeline_lake_fungi/rename.pl
cat renamed*.fas > "${TMP_FASTA1}"
rm renamed*.fas

# Dereplicate (vsearch); removing the option --relabel_sha1 \
"${VSEARCH}" --derep_fulllength "${TMP_FASTA7}" \
             --sizein \
             --sizeout \
             --minseqlength 10 \
             --fasta_width 0 \
	     --minuniquesize ${LOW_ABUND_SAMPLE} \
	     --sizein -sizeout \
             --threads ${THREADS} \
             --output "${TMP_FASTA2}" > /dev/null

# Sorting
"${VSEARCH}" --fasta_width 0 \
             --sortbysize "${TMP_FASTA2}" \
             --minseqlength 10 \
             --threads ${THREADS} \
	--output $PROJ.centroids --sizein --sizeout


rm -f $PROJ.uc "${TMP_FASTA1}" "${TMP_FASTQ2}" "${TMP_FASTQ7}"