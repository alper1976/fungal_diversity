######################################
### PIPELINE LAKE FUNGI 
###
###  Pedro Martin-Sanchez - July 2022
#########################

######### RELEVANT DATA of the 2 projects
 
####### 1) NORWAY
# 
# 74 lake samples, failure is expected in one of them (L24; no band in the second PCR)
# LNegPCR with evident contamination in the second PCR
# LMock: 3 fungal species, used in previous projects  
#
# PRIMERS USED (1XF and 2xR- equimolar):
# ITS3-Mix2 (forward)
# 5’ -ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNN + caWCGATGAAGAACGCAg -3’
#
# ITS4-cwmix1 (reverse)
# 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNN + TCCTCCGCTTAyTgATAtGc -3’ 
# ITS4-cwmix2 (reverse)
# 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNN + TCCTCCGCTTAtTrATAtGc -3’
#
# So look for the adapters:
# FWD: CAWCGATGAAGAACGCAG
# REV: TCCTCCGCTTAYTRATATGC (2 wildcard characters: Y and R)
#
# But also for:  
# FWD.RC:  CTGCGTTCTTCATCGWTG
# REV.RC: GCATATYARTAAGCGGAGGA

####### 2) SWEDEN

See run file (run021.xlsx)
See metadata file (100lakes_metadata_final.xlsx)

After a first look, this seems to contain 177 samples:		##Keep all samples for the bioinformatics
- 103 samples from Swedish lakes
- 28 samples from Canada
- 22+24 samples from Swedish and Finnish lakes with variations (different zones, pore diameters...?)
 
PRIMERS (Extracted from run021.xlsx):		# Almost identical primers, so identical target to the norwegian study

ITS3mix2	AACCAWCGATGAAGAACGCAG		#Extra AAC in 5' compared to the other study (CAWCGATGAAGAACGCAG)
ITS4cwmix	AATCCTCCGCTTAYTRATATGC		#Extra AA in 5' compared to the other study (TCCTCCGCTTAYTRATATGC)

So look for the adapters:
FWD: AACCAWCGATGAAGAACGCAG
REV: AATCCTCCGCTTAYTRATATGC (2 wildcard characters: Y and R)
     
But also for:  
FWD.RC: CTGCGTTCTTCATCGWTGGTT
REV.RC: GCATATYARTAAGCGGAGGATT

####### BOTH PROJECT:  TARGET AMPLICON: about 350-500pb including aprox. 130bp 5.8S rRNA gene and the full ITS2 region (highly variable in both sequence and length)   

# Similar target was used by : Heeger F, Wurzbacher C, Bourne EC, Mazzoni CJ, Monaghan MT. 
# Combining the 5.8S and ITS2 to improve classification of fungi. Methods Ecol Evol. 2019;10:1702â€“1711. 
# https://doi.org/10.1111/2041-210X.1326
#
# NB: They developed an specific pipeline using a python snake - see github
# https://github.com/f-heeger/two_marker_metabarcoding

# I tried to adapt it to a similar pipeline using DADA2 among other changes... 


##########################################################
### PIPELINE SUMMARY
##########################

########## INITIAL PREPARATIONS/INSTALLATIONS:
#### Install the needed R packages including:
# devtools  
install.packages("devtools", repos="https://cran.uib.no/")
# dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10") 
# lulu
install_github("tobiasgf/lulu")
#openxlsx to save some tables from the DADA2 script
install.packages("openxlsx", dependencies=TRUE, repos="https://cran.uib.no/")

#### create a conda environment including  other software needed in some steps of the pipeline: cutadapt, itsx, fasplit, vsearch
module purge
module load Anaconda3/2019.03
conda create -n fungi_metabarcoding_pipeline -c bioconda -c conda-forge cutadapt=2.7 itsx ucsc-fasplit vsearch=2.14.0 --yes 
conda init bash

########## The developed pipeline includes the following steps
 (the script mentioned are those used for the project 100lakes_Sweden):

######## 1) CUTADAPT: trimming the primers (also looking for reverse complementary ones) and removal of reads <200 bp

# Run interactively in Saga
script_cutadapt_sweden


######## 2) DADA2, with some changes for fungi (checking presence of primers and NO truncLen), includes: (i) Quality filtering and trimming (Phred quality>20 ->truncQ=20; minLen=50;max_n=0; maxEE=2),
   (ii) Dereplication, denoising and merging, and(iii) Removal of chimeras (42,7% of ASVs!!) -> 2458 ASVs
# Adaptation of DADA2 outputs (create a folder with fasta files for each sample, edit samples names, etc.) for the next step...

# Run the slurm file: 
slurm_DADA2

# associated with the R code:
sweden_lake_fungi_dada2_Rcode.R
  
# Generated OUTPUTS:
# folder with extracted reads in fasta files (one per sample): "DADA2_extracted_samples_no_chim" 
# "DADA2_noChime.table" where the rows are ASVs (Seq1-Seq226) and columns are samples with namees like: 25-L41_S25_L001.fastq.gz) 
# "DADA2_noChime.raw.otus": sequences in fasta format, like: >Seq1.... 


######## 3) ITSx - extraction of regions ITS2 and 5.8S for ALL organisms (No filter fungi)
# Activate the Conda environment with the needed softwares
conda activate fungi_metabarcoding_pipeline

# Run interactively from the DADA2_extracted_samples_no_chim the corresponding script:
bash script_parallel_itsx_fungi

# This created 2 folders with trimmed sequences in fasta files for each sample. "ITS2"  and "5_8S"


######## 4.1) for ITS2:

######## 4.1.1 - Clustering at 97% and removal of singletons using VSEARCH 

# Activate the Conda environment with the needed softwares
conda activate fungi_metabarcoding_pipeline

# Run interactively the corresponding code from the folder DADA2_extracted_samples_no_chim/ITS2
bash script_VSEARCH_cluster_sweden.sh
# The script should be adapted accordingly depending on sample names!!
# The exchange of strings using sed, considering the hetereogeneous sample labels, could be problematic intruducing errors - we should double check later!  

# Output files:
# Aqua_Fungi_ITS2_Sweden.otutable 
# Aqua_Fungi_ITS2_Sweden.centroids (the representative ITS2 sequences)
# Aqua_Fungi_ITS2_Sweden.uc (I kept this intermediate file to track the changes done in the preparation for making the otutable,

# The samples names in the otutable files should be corrected manualy!!
# Delete "barcodelabel=" from all sample names (I used the tool replace in the WinSCP editor)  


######## 4.1.2 - Curation by LULU -> final OTU table

# Run the slurm file from the folder DADA2_extracted_samples_no_chim/ITS2:
bash script_6_LULU_corrected.sh

# Associated with the corresponding R code:
scriptR_code_LULU

# This was very quick and created the following output files associated with database and LULU results:
# OTU_centroids2 (7 files)
# match_list.txt
# lulu.log_xxxxxxx
# lulu.analyses.R (rds file)
# OTU_table_lulu_curated (the final curated otutable) 

## Before BLAST and LCA, we shift to 5.8S....

######## 4.2) for 5.8S:

######## 4.2.1) - Dereplication, fitering (--minseqlength 10, --minuniquesize 2) and sorting with VSEARCH
          (Input: fasta files for each sample after DADA2 and ITSx -> Output: an unique fasta file with 1296 ASV representative sequences)    

# To get the needed "centroids" file for BLAST
# In the folder "DADA2_extracted_samples_no_chim/5_8S"

# Activate the Conda environment with the needed softwares
conda activate fungi_metabarcoding_pipeline

# Run interactively the corresponding code:
bash script_VSEARCH_derep_sort_filtering_Sweden.sh

# I got the needed output: Aqua_Fungi_5_8S.centroids
# Note: they are not really centroids! as there was no clustering here!


########  TAXONOMIC ASSINGMENT by BLAST, separately for ITS2 and 5.8S 
# From the folders "DADA2_extracted_samples_no_chim/ITS2" or "DADA2_extracted_samples_no_chim/5_8S"

######### 4.1.3) ITS2 database (in ITS2 folder)

module load BLAST+/2.8.1-foss-2018b

# (Already done) prepare the database from the fasta file (ITS2)
#makeblastdb -in /cluster/projects/nn9745k/03_databases/fungi/sh_general_release_all_10.05.2021/ITS2_unite_insd_2021_all.fas -out /cluster/projects/nn9745k/03_databases/fungi/sh_general_release_all_10.05.2021/ITS2_unite_insd_2021_all_db -dbtype nucl -title "UNITE ITS2 database"
## "Adding sequences from FASTA; added 101690 sequences in 2.80823 seconds." NOTE: 9000 SEQUENCES LESS after ITSx!!

mkdir blast_lca_test

# Run BLAST
blastn -max_target_seqs 100 -evalue 1 -query /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/ITS2/Aqua_Fungi_ITS2_Sweden.centroids -out /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/ITS2/blast_lca_test/blast_results_ITS2 -db /cluster/projects/nn9745k/03_databases/fungi/sh_general_release_all_10.05.2021/ITS2_unite_insd_2021_all_db -outfmt 6 -num_threads 2


######### 4.2.2) 5.8S database (in 5_8S folder)

#(Already done) prepare the database from the fasta file (5.8S)
#makeblastdb -in /cluster/projects/nn9745k/03_databases/fungi/sh_general_release_all_10.05.2021/5.8S_unite_insd_2021_all.fas -out /cluster/projects/nn9745k/03_databases/fungi/sh_general_release_all_10.05.2021/5.8S_unite_insd_2021_all_db -dbtype nucl -title "UNITE 5.8S database" 
### "Adding sequences from FASTA; added 101394 sequences in 3.16221 seconds.." NOTE: about 9000 SEQUENCES LESS after ITSx!!

mkdir blast_lca_test

# Run BLAST
blastn -max_target_seqs 100 -evalue 1 -query /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/5_8S/Aqua_Fungi_5_8S_Sweden.centroids -out /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/5_8S/blast_lca_test/blast_results_5_8S -db /cluster/projects/nn9745k/03_databases/fungi/sh_general_release_all_10.05.2021/5.8S_unite_insd_2021_all_db -outfmt 6 -num_threads 2


###### We got the 2 needed outputs for the next steps:
#blast_results_ITS2  
#blast_results_5_8S


####### LCA analyses, separately for ITS2 and 5.8S 

# first we need to adapt the formating of the blast results for subsequent LCA anayses
# using the following R code in Saga:
Adapt_blast_lca.R


######## LCA analyses in Bash: 

# Move the lca.py script to the folder: /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi
# Selected PARAMETERS: -b 8 -id 80 -cov 85 -t only_lca -flh unidentified

######## 4.1.4) ITS2
python lca.py -i /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/ITS2/blast_lca_test/blast_results_for_lca_ITS2 -o /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/ITS2/blast_lca_test/blast_results_lca_8_80_85_only_lca_flh_ITS2 -b 8 -id 80 -cov 85 -t only_lca -flh unidentified

######## 4.2.3) 5.8S
python lca.py -i /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/5_8S/blast_lca_test/blast_results_for_lca_5_8S -o /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi/DADA2_extracted_samples_no_chim/5_8S/blast_lca_test/blast_results_lca_8_80_85_only_lca_flh_5_8S -b 8 -id 80 -cov 85 -t only_lca -flh unidentified


######## 5) Final taxonomic classification (combination of ITS2 <- 5.8S) 

#### As described by Heeger et al 2019, the purpose in this final classification is:

"2.3-d) The final classification of the OTUs defined by the ITS2 combines the classifications from stage 2 and 3 (#LCAs fron ITS2 and 5.8S, I guess). 
For each read present in an ITS2 OTU cluster, all 5.8S sequences and their classifications are collected.
The 5.8S classifications are combined with the same LCA approach explained above. 
The resulting classification is compared to the ITS2 classification.
If 5.8S and ITS2 classification are concordant, but the ITS2 is classified to a lower taxonomic rank, the ITS2 classification is accepted.
Sequences that are unclassified with ITS2 will automatically take the 5.8S classifications. 
All conflicting classifications can either be marked (default) or resolved by the user by giving priority to one of the marker." 

### The selected files for the final classification were moved to the folder "final_tax_classification":
- Final LULU-curated OTU table including all organisms: "OTU_table_lulu_curated.xlsx"
- the selected LCA results for ITS2 
- the selected LCA results for 5.8S

# Run R script for merging two taxonomies in Saga:"merge_taxonomy.R"  

### This saves the following output files:
# LCA merged taxonomy
# final OTU table with taxonomy (xlsx) 
# phyloseq objects including OTUs and taxonomy

################## END ####################