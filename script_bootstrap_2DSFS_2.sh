#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
#SBATCH --array=1-100  ###3021Main chromosomes only in the intervals_output.txt file
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=8
# Walltime format hh:mm:ss
#SBATCH --time=03:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module load angsd/0.940
#realSFS Results_for_load/countryside.list.saf.idx  Results_for_load/urban_unrelated.list.saf.idx -sites LOW_SNPs.sites -bootstrap 1 -P 10 -m 0 > Results_for_load/country_urban_LOW_boot_${SLURM_ARRAY_TASK_ID}.ml


realSFS Results_for_load/countryside.list.saf.idx  Results_for_load/urban_unrelated.list.saf.idx -sites MEDIUM_SNPs.sites -bootstrap 1 -resample_chr 1 -seed ${SLURM_ARRAY_TASK_ID} -P 10 -m 0 > Results_for_load/country_urban_MEDIUM_boot_${SLURM_ARRAY_TASK_ID}.ml


#realSFS Results_for_load/countryside.list.saf.idx  Results_for_load/urban_unrelated.list.saf.idx -sites HIGH_SNPs.sites -bootstrap 1 -P 10 -m 0 > Results_for_load/country_urban_HIGH_boot_${SLURM_ARRAY_TASK_ID}.ml
