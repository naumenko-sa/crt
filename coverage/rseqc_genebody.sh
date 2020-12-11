#!/bin/bash

# calculate gene body coverage to see if it is skewed towards 3UTR
# http://rseqc.sourceforge.net/#genebody-coverage-py

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=short          # Partition (queue)
#SBATCH --time=12:00:00            # Runtime in D-HH:MM format
#SBATCH -c 1
#SBATCH --job-name=rseqc            # Job name
#SBATCH --mem=10G                   # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

module load gcc/6.2.0
module load python/2.7.12
module load rseqc/2.6.4

# $1 = exons.bed
# $2 = bam.list

geneBody_coverage.py -r $1 -i $2 -o seqc_output

