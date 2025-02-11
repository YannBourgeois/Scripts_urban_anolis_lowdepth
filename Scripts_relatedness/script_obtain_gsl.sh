#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
#SBATCH --array=1-1000  ###Main chromosomes only in the intervals_output.txt file
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Walltime format hh:mm:ss
#SBATCH --time=24:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module load angsd/0.940

config=intervals_output.txt
REGION=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)


#angsd -GL 1 -out relatedness_from_BAM/genolike_${SLURM_ARRAY_TASK_ID} -nThreads 10 -r ${REGION} -doGlf 3 -ref anoCar2.fa -doMajorMinor 1 -minmaf 0.05 -SNP_pval 1e-8 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 55 -bam list_all_New_Orleans.txt
#angsd -b list_all_demo.list -gl 2 -domajorminor 1 -rf chromo.list -snp_pval 1e-8 -domaf 1 -minmaf 0.05 -doGlf 3 -nThreads 50 -out genolike
angsd -GL 1 -out relatedness_from_BAM/genolike_${SLURM_ARRAY_TASK_ID} -nThreads 10 -r ${REGION} -doGlf 2 -ref anoCar2.fa -doMajorMinor 1 -minmaf 0.05 -SNP_pval 1e-8 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 55 -bam list_all_New_Orleans.txt
