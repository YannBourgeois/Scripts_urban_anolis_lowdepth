#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
#SBATCH --array=1-110
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

module load samtools
module load bedtools


config=array_bedtools.txt
BAM=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

bedtools genomecov -bga -ibam ${BAM} -g scaffolds_length.txt > ${BAM}_coverage.bed
awk -F"\t" '$4==0 || $4>10' ${BAM}_coverage.bed > ${BAM}_mask_10x.bed
awk -F"\t" '$4==0 || $4>20' ${BAM}_mask_10x.bed > ${BAM}_mask_20x.bed
awk -F"\t" '$4==0 || $4>25' ${BAM}_mask_20x.bed > ${BAM}_mask_25x.bed
