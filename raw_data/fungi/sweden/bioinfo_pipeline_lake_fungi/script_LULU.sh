module purge
module load BLAST+/2.8.1-foss-2018b
OTU_centroids=$(find -type f -name '*.centroids')
sed 's/;.*$//' "${OTU_centroids}" > OTU_centroids2
makeblastdb -in OTU_centroids2 -parse_seqids -dbtype nucl
blastn -db OTU_centroids2 -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query OTU_centroids2

# run the LULU curation
module purge
module load R/3.5.1-intel-2018b
R < script_R_code_LULU --no-save

