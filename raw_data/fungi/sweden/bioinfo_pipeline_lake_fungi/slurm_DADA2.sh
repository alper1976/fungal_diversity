#!/bin/bash

#SBATCH --account=nn9744k
#SBATCH --job-name=dada2
#SBATCH --time=1-12:0:0
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=12G


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module purge
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
R < sweden_lake_fungi_dada2_Rcode.R --no-save

