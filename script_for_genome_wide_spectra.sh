#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
#SBATCH --array=1-1007  ###3021Main chromosomes only in the intervals_output.txt file
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=10
# Walltime format hh:mm:ss
#SBATCH --time=02:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module load angsd/0.940

realSFS Results_for_FST_final/countryside_intervals_${SLURM_ARRAY_TASK_ID}.saf.idx Results_for_FST_final/urban_unrelated_intervals_${SLURM_ARRAY_TASK_ID}.saf.idx -P 10 -m 0 -sites INTERGENIC_SNPs.sites INTERGENIC_SNPs.sites > Results_for_load/country_urban_GENOMEWIDE_interval_${SLURM_ARRAY_TASK_ID}.ml
