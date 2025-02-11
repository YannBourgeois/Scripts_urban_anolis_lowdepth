#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
#SBATCH --array=1-71  ###3021Main chromosomes only in the intervals_output.txt file
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
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
config=array_SFS_individual_BAM.txt
BAM=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

angsd -i ${BAM} -doCounts 1 -doMaf 1 -doMajorMinor 3 -doSaf 1 -GL 1 \
-sites INTERGENIC_SNPs.sites -nThreads 10 -ref anoCar2.fa  \
-anc anoCar2.fa  \
-out Results_for_individual_load/${BAM}_INTERGENIC -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 30


realSFS Results_for_individual_load/${BAM}_INTERGENIC.saf.idx -sites INTERGENIC_SNPs.sites -P 4 -m 0 > Results_for_individual_load/${BAM}_INTERGENIC.ml

