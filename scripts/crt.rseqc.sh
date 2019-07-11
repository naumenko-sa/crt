#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority          # Partition (queue)
#SBATCH --time=1-00:00              # Runtime in D-HH:MM format
#SBATCH -c 1
#SBATCH --job-name=sahay            # Job name
#SBATCH --mem-per-cpu=10G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

module load gcc/6.2.0
module load python/2.7.12
module load rseqc/2.6.4

geneBody_coverage.py -r GRCm38_mm10_Ensembl.nochr.bed -i gene_body -o output

