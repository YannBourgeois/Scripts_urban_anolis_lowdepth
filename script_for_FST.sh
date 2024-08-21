#!/bin/bash
#SBATCH --mem=16GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
#SBATCH --array=1-1007 ###3021Main chromosomes only in the intervals_output.txt file
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

config=array3.txt
SUBSET0=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
SUBSET1=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
SUBSET2=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

./winsfs Results_for_FST_final/${SUBSET0}.saf.idx  Results_for_FST_final/${SUBSET1}.saf.idx > Results_for_FST_final/${SUBSET0}_${SUBSET2}.ml
./winsfs Results_for_FST_final/${SUBSET0}.saf.idx  Results_for_FST_final/${SUBSET2}.saf.idx > Results_for_FST_final/${SUBSET0}_${SUBSET1}.ml
./winsfs Results_for_FST_final/${SUBSET1}.saf.idx  Results_for_FST_final/${SUBSET2}.saf.idx > Results_for_FST_final/${SUBSET1}_${SUBSET2}.ml

grep -v "#" Results_for_FST_final/${SUBSET0}_${SUBSET2}.ml > Results_for_FST_final/${SUBSET0}_${SUBSET2}.ml2;mv Results_for_FST_final/${SUBSET0}_${SUBSET2}.ml2 Results_for_FST_final/${SUBSET0}_${SUBSET2}.ml
grep -v "#" Results_for_FST_final/${SUBSET0}_${SUBSET1}.ml > Results_for_FST_final/${SUBSET0}_${SUBSET1}.ml2;mv Results_for_FST_final/${SUBSET0}_${SUBSET1}.ml2 Results_for_FST_final/${SUBSET0}_${SUBSET1}.ml
grep -v "#" Results_for_FST_final/${SUBSET1}_${SUBSET2}.ml > Results_for_FST_final/${SUBSET1}_${SUBSET2}.ml2;mv Results_for_FST_final/${SUBSET1}_${SUBSET2}.ml2 Results_for_FST_final/${SUBSET1}_${SUBSET2}.ml

realSFS fst index Results_for_FST_final/${SUBSET0}.saf.idx  Results_for_FST_final/${SUBSET1}.saf.idx Results_for_FST_final/${SUBSET2}.saf.idx -sfs Results_for_FST_final/${SUBSET0}_${SUBSET2}.ml -sfs Results_for_FST_final/${SUBSET0}_${SUBSET1}.ml -sfs Results_for_FST_final/${SUBSET1}_${SUBSET2}.ml -fstout Results_for_FST_final/FST_${SLURM_ARRAY_TASK_ID}_fst -whichFst 1
realSFS fst stats2 Results_for_FST_final/FST_${SLURM_ARRAY_TASK_ID}_fst.fst.idx -win 50000 -step 10000 -type 0 >  Results_for_FST_final/FST_${SLURM_ARRAY_TASK_ID}_50kb.fst2
realSFS fst stats2 Results_for_FST_final/FST_${SLURM_ARRAY_TASK_ID}_fst.fst.idx -win 5000 -step 1000 -type 0 >  Results_for_FST_final/FST_${SLURM_ARRAY_TASK_ID}_5kb.fst2

