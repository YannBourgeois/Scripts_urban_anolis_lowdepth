#!/bin/bash
#SBATCH --mem=10GB
#SBATCH --job-name=RSYNC
#SBATCH --time=72:00:00
#SBATCH --partition long
#SBATCH --nodes=1
#SBATCH --output=rsync.%J.out
#SBATCH --error=rsync.%J.err
mv Coffea_project/ /shared/projects/datepalm/

