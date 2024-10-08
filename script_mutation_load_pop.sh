#!/bin/bash
#SBATCH --mem=64GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
#SBATCH --array=2-3  ###3021Main chromosomes only in the intervals_output.txt file
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=8
# Walltime format hh:mm:ss
#SBATCH --time=36:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module load angsd/0.940

config=array_2.txt

POP=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
temp=`wc -l ${POP} | cut -f 1 -d " "`
NIND=$((2*${temp}/3)) ###at least 66% individuals sequenced

if [ "${POP}" = "GA.list" ]; then
    maxDepth=25
elif [ "${POP}" = "countryside.list" ]; then
    maxDepth=20
elif [ "${POP}" = "urban_unrelated.list" ]; then
    maxDepth=10
else
    maxDepth=0  # Default value if no conditions are met
fi

angsd -doCounts 1 -doSaf 2 -sites Annotated_SNPs.sites -GL 1 \
-nThreads 8 -ref anoCar2.fa  \
-anc anoCar2.fa  \
-doMaf 1 -doMajorMinor 3 -minInd ${NIND} -setMaxDepth $(($temp*$maxDepth)) -rf main_chromosomes.txt -out Results_for_load/${POP} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 30 -bam ${POP} -indF ${POP}.indF

#realSFS Results_for_load/${POP}.saf.idx -P 8 -sites LOW_SNPs.sites -m 0 > Results_for_load/${POP}.ml
#realSFS Results_for_load/${POP}.saf.idx -P 8 -sites MEDIUM_SNPs.sites -m 0 > Results_for_load/${POP}.ml
#realSFS Results_for_load/${POP}.saf.idx -P 8 -sites HIGH_SNPs.sites -m 0 > Results_for_load/${POP}.ml
