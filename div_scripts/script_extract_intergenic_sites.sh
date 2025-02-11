#!/bin/bash
#SBATCH --mem=32GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Walltime format hh:mm:ss
#SBATCH --time=12:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
rm -r ../datepalm/Results_for_FST_final/
rm -r ../datepalm/Results_for_FST_without_inbreeding_filters/
mv Results_for_FST_final/ ../datepalm/
grep intergenic_region all_filtered_ann.vcf  | cut -f1,2,4,5 | grep -v "," | sort -k1,1 -k2,2n > INTERGENIC_SNPs.sites
grep intergenic_region small_filtered_ann.vcf  | cut -f1,2,4,5 | grep -v "," | sort -k1,1 -k2,2n >> INTERGENIC_SNPs.sites

shuf -n 100000 INTERGENIC_SNPs.sites | sort -k1 -k2,2n > SUBSET_INTERGENIC_SNPs.sites
