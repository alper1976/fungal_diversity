#!/bin/bash
# Job name: tree fungal
#SBATCH --job-name=fungi_tree
#
# Project (change to your :
#SBATCH --account=nn9745k
#
# Wall clock limit:
#SBATCH --time=10:0:0
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=8000MB

#
## Number of tasks (nodes):
#SBATCH --ntasks-per-node=8
##


##### setting up environment

module load MAFFT/7.490-GCC-11.2.0-with-extensions
module load FastTree/2.1.11-GCCcore-11.2.0

# This scripts performs sequence alignments and makes phylogenetic trees for community assembly analyses

# Norwegian dataset

cd /cluster/projects/nn9744k/02_results/29_euk_20220325/Pedro

mafft --maxiterate 1000 --thread 8 --globalpair rep-seqs.fna > rep-seqs_aligned.fna
FastTree -nt -nosupport -fastest rep-seqs_aligned.fna > rep-seqs.tree

# Swedish dataset

cd /cluster/projects/nn9745k/02_results/02_hundred_lakes/fungi

mafft --maxiterate 1000 --thread 8 --globalpair rep-seqs.fna > rep-seqs_aligned.fna
FastTree -nt -nosupport -fastest rep-seqs_aligned.fna > rep-seqs.tree