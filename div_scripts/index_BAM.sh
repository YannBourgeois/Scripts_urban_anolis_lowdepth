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
#module load bcftools/1.17 


config=list_to_index.txt
BAM=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

samtools index ${BAM}
