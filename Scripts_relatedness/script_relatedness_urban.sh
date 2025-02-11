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

#vcftools --vcf genptyped_yann_cgvcf.vcf --max-missing 0.6 --maxDP 10 --minQ 30 --minGQ 10 --keep to_keep_demography_vcf.txt --mac 8 --recode --recode-INFO-all --out subset_for_ngrRelate 
#./BAM_files/ngsRelate/ngsRelate -h subset_for_ngrRelate.recode.vcf -O relatedness.res -T GT -c 1 -p 20 -z individuals_demo_ID.txt -n 80
#vcftools --vcf genptyped_yann_cgvcf.vcf --maxDP 10 --minQ 30 --keep to_keep_demography_vcf.txt --mac 8 --recode --recode-INFO-all --out subset_for_ngrRelate

##../ngsRelate/ngsRelate -h all_filtered.vcf -O relatedness_VCF.res -p 30 -n 103

vcftools --vcf all_filtered.vcf --maxDP 10 --minQ 30 --mac 8 --keep low_missingness_urban.txt --recode --recode-INFO-all --out ngsRelate_urban
./ngsRelate/ngsRelate -h ngsRelate_urban.recode.vcf -O relatedness_VCF_urban.res -p 30 -z low_missingness_urban.txt
