#!/bin/bash
#SBATCH --mem=64GB
#SBATCH -p fast
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=5
# Walltime format hh:mm:ss
#SBATCH --time=06:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
#module load gcc/4.9.3
#module load gencore/1
#module load gencore_variant_detection/1.0
#module load htslib/

module load vcftools/

vcftools --vcf all_filtered.vcf --maxDP 10 --minQ 30 --mac 1 --keep low_missingness_country.txt --recode --recode-INFO-all --out ngsRelate_country
./ngsRelate/ngsRelate -h ngsRelate_country.recode.vcf -O relatedness_VCF_country.res -p 30 -z low_missingness_country.txt

