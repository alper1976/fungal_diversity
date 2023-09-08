module purge

#SETTING VARIABLES
VSEARCH=$(which vsearch)
THREADS=1
CLUSTER=0.97
LOW_ABUND_SAMPLE=2
PROJ="Aqua_Fungi_ITS2_Sweden"
TMP_FASTA1=$(mktemp)
TMP_FASTA2=$(mktemp)
TMP_FASTA3=$(mktemp)
TMP_FASTA4=$(mktemp)
TMP_FASTA5=$(mktemp)
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
             --output "${TMP_FASTA3}" --sizein --sizeout

#Clustering
"${VSEARCH}" --cluster_size "${TMP_FASTA3}" \
	     --id ${CLUSTER} \
	     --sizein --sizeout \
             --minseqlength 10 \
             --qmask none \
             --threads ${THREADS} \
	     --centroids "${TMP_FASTA4}"

# Sorting & filtering
"${VSEARCH}" --fasta_width 0 \
             --minseqlength 10 \
             --sortbysize "${TMP_FASTA4}" \
             --threads ${THREADS} \
             --output $PROJ.centroids --sizein --sizeout

#mapping of raw reads against OTU representatives
"${VSEARCH}" --usearch_global "${TMP_FASTA1}" \
	     --db $PROJ.centroids  \
             --minseqlength 10 \
	     --strand plus \
	     --id ${CLUSTER} \
	     --maxaccepts 0 \
             --qmask none \
             --threads ${THREADS} \
	     --uc $PROJ.uc

#preparing for making table
sed -i 's/	A/	barcodelabel=A/g' $PROJ.uc
sed -i 's/	B/	barcodelabel=B/g' $PROJ.uc
sed -i 's/	C/	barcodelabel=C/g' $PROJ.uc
sed -i 's/	D/	barcodelabel=D/g' $PROJ.uc
sed -i 's/	E/	barcodelabel=E/g' $PROJ.uc
sed -i 's/	F/	barcodelabel=F/g' $PROJ.uc
sed -i 's/	g/	barcodelabel=g/g' $PROJ.uc
sed -i 's/	G/	barcodelabel=G/g' $PROJ.uc
sed -i 's/	H/	barcodelabel=H/g' $PROJ.uc
sed -i 's/	I/	barcodelabel=I/g' $PROJ.uc
sed -i 's/	L/	barcodelabel=L/g' $PROJ.uc

sed -i 's/	1/	barcodelabel=1/g' $PROJ.uc
sed -i 's/	2/	barcodelabel=2/g' $PROJ.uc
sed -i 's/	3/	barcodelabel=3/g' $PROJ.uc
sed -i 's/	4/	barcodelabel=4/g' $PROJ.uc
sed -i 's/	5/	barcodelabel=5/g' $PROJ.uc
sed -i 's/	6/	barcodelabel=6/g' $PROJ.uc
sed -i 's/	7/	barcodelabel=7/g' $PROJ.uc
sed -i 's/	8/	barcodelabel=8/g' $PROJ.uc
sed -i 's/	9/	barcodelabel=9/g' $PROJ.uc


#make table
python2 /cluster/projects/nn9745k/scripts/bioinfo_pipeline_lake_fungi/uc2otutab.py $PROJ.uc > $PROJ.otutable

rm -f "${TMP_FASTA1}" "${TMP_FASTQ2}" "${TMP_FASTQ3}" "${TMP_FASTQ4}" "${TMP_FASTQ5}" "${TMP_FASTQ7}"
