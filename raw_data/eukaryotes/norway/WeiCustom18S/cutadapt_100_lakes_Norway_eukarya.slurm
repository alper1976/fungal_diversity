#!/bin/bash
# Job name: eukaryota_cutadapt
#SBATCH --job-name=e_cutadapt
#
# Project (change to your :
#SBATCH --account=nn9745k
#
# Wall clock limit:
#SBATCH --time=2:0:0
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=32G
#
## Number of tasks (nodes):
#SBATCH --ntasks=1


##Setup job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

module load cutadapt/1.18-foss-2018b-Python-3.6.6

######################################
### Eukarya hundred_lakes project

## Illumina run WeiCustom18S

# mkdir /cluster/projects/nn9745k/02_results/03_hundred_lakes_norway
# mkdir /cluster/projects/nn9745k/02_results/03_hundred_lakes_norway/eukaryota
# mkdir /cluster/projects/nn9745k/02_results/03_hundred_lakes_norway/eukaryota/AdaptersRemoved

##############################################
### run cutadapt to remove primer sequences
##############################################

DATADIR=/cluster/projects/nn9745k      #Root of all data
INPUT=$DATADIR/01_raw_data/03_hundred_lakes_norway/WeiCustom18S
OUTPUT=$DATADIR/02_results/03_hundred_lakes_norway/eukaryota/AdaptersRemoved

for f1 in $INPUT/*_R1_*.fastq.gz
    do
        f2=${f1/_R1_/_R2_}
        out1=$OUTPUT/$(basename ${f1})
        out2=${out1/_R1_/_R2_}
        cutadapt -g CCAGCASCYGCGGTAATTCC -G ACTTTCGTTCTTGATYRATGA -m 200 -e 0 --discard-untrimmed -o ${out1} -p ${out2} ${f1} ${f2}
    done

## End of script
