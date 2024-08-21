#!/bin/bash
#SBATCH --mem=64GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=32
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

conda activate ngsF
LD_LIBRARY_PATH=/shared/home/ybourgeois/.conda/envs/ngsF/lib/
export LD_LIBRARY_PATH

NSAMS=`cat list_all_New_Orleans.txt | wc -l`
NSITES=`cat relatedness_from_BAM/freq_partial_ngsrelate | wc -l`

#./ngsF/ngsF.sh --n_ind $NSAMS --n_sites $NSITES --n_threads 10 --glf relatedness_from_BAM/combined_genolike.glf --out Inbreeding_ngsF.indF

./ngsF-HMM/ngsF-HMM.sh --n_ind $NSAMS --n_sites $NSITES --n_threads 32 --geno relatedness_from_BAM/combined_genolike_ngsFHMM.glf --pos relatedness_from_BAM/positions_ngsF-HMM.txt --out Inbreeding_ngsF_HMM
