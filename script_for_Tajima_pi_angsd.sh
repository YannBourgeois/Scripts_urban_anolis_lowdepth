#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
#SBATCH --array=1-3021  ###3021Main chromosomes only in the intervals_output.txt file
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
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

config=array.txt
SUBSET=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
REGION=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
POP=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
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

angsd -doCounts 1 -doSaf 2 -doSnpStat 1 -doHWE 1 -rf ${REGION} -GL 1 \
-nThreads 10 -ref anoCar2.fa  \
-anc anoCar2.fa  \
-doMaf 1 -doMajorMinor 1 -minInd ${NIND} -out Results_for_FST_final/${SUBSET} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 30 -setMaxDepth $(($temp*$maxDepth)) -bam ${POP} -indF ${POP}.indF
realSFS Results_for_FST_final/${SUBSET}.saf.idx -P 4 -nSites 2000000 -m 0 > Results_for_FST_final/${SUBSET}.ml
realSFS saf2theta Results_for_FST_final/${SUBSET}.saf.idx -sfs Results_for_FST_final/${SUBSET}.ml -outname Results_for_FST_final/${SUBSET}
thetaStat do_stat Results_for_FST_final/${SUBSET}.thetas.idx -win 50000 -step 10000 -outnames Results_for_FST_final/${SUBSET}_50kb.thetas
thetaStat do_stat Results_for_FST_final/${SUBSET}.thetas.idx -win 5000 -step 1000 -outnames Results_for_FST_final/${SUBSET}_5kb.thetas
