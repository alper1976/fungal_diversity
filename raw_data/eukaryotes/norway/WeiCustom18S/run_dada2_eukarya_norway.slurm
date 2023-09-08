#!/bin/bash
# Job name: dada2 Eukarya Norway 100 lakes
#SBATCH --job-name=dada2_e
#
# Project (change to your :
#SBATCH --account=nn9745k
#
# Wall clock limit:
#SBATCH --time=20:0:0
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=61G

#
## Number of tasks (nodes):
##SBATCH --ntasks-per-node=12
##

##Setup job environment
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

#Load R module

module load R/3.6.0-foss-2019a

#Set variable with number of processors and load mpi

mpirun -np 1 R CMD BATCH --no-save --no-restore /cluster/projects/nn9745k/scripts/03_hundred_lakes_norway/eukarya/dada2_script_100_lakes_Norway_eukarya.R

#End of script
