#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
####SBATCH --array=22-3021  ###3021Main chromosomes only in the intervals_output.txt file
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=10
# Walltime format hh:mm:ss
#SBATCH --time=12:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module load angsd/0.940

SUBSET=GA
#realSFS cat -b SAF_GA.list -outnames GA -P 10
realSFS ${SUBSET}.saf.idx -P 4 -sites LOW_SNPs.sites -m 0 > Results_for_load_2/${SUBSET}_LOW.ml
realSFS ${SUBSET}.saf.idx -P 4 -sites MEDIUM_SNPs.sites -m 0 > Results_for_load_2/${SUBSET}_MEDIUM.ml
realSFS ${SUBSET}.saf.idx -P 4 -sites HIGH_SNPs.sites -m 0 > Results_for_load_2/${SUBSET}_HIGH.ml
