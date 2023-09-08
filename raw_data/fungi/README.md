# raw_data
These are scripts (slurm scripts and R code) used to run raw data processing on sequencing runs 131002_M00485_0074_000000000-A5GKR, 160408_M00485_0264_000000000-AME7C and 161010_M00485_0309_000000000-ATHUN. 

## individual scripts
### runscript_dada2.sh
This contains bash code to demultiplex 131002_M00485_0074_000000000-A5GKR and run cutadapt on all sequencing runs.

### run_dada2_abel.slurm 
This is a slurm script that runs the R script dada2_abel.R

### dada2_abel.R 
This is an R script that runs dada2 to obtain ASV tables and taxonomic annotations.

