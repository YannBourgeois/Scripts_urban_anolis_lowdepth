#!/bin/bash
#SBATCH --mem=32GB
#SBATCH --job-name=arg
#SBATCH --time=72:00:00
#SBATCH --partition long
#SBATCH --array=1-40
#SBATCH --nodes=1
#SBATCH --output=arg.%J.out
#SBATCH --error=arg.%J.err

module load htslib
module load vcftools

#vcftools --vcf all_filtered.vcf --keep list_to_keep_29unrelated_depth2x.txt --mac 1 --recode --recode-INFO-all --chr chr4 --out chr4_filtered

#bgzip -c chr4_filtered.recode.vcf > chr4_filtered.recode.vcf.gz

#tabix -p vcf chr4_filtered.recode.vcf.gz
config=array_arg.txt
REGION=$(awk -F"\t" -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

conda activate ARG

arg-sample --vcf chr4_filtered.recode.vcf.gz \
--popsize-file Population_sizes_rescaled_higher_rate.txt \
--ind-maskmap list_indv_masks.txt \
--unphased \
--mask-Ns 52 --mask-cluster 2,5 --vcf-min-qual 30 -t 15 \
--use-genotype-probs \
-r 6.68e-10 -m 1.7e-9 \
--maxtime 5.00E+05 \
-n 2000 \
--region ${REGION} -o chr4_${SLURM_ARRAY_TASK_ID}

