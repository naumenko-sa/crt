#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=kallisto          # Job name
#SBATCH -c 8			    # cores
#SBATCH --mem=10G           # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date
which kallisto

# https://bustools.github.io/BUS_notebooks_R/velocity.html

# 100G ram
#kallisto bus \
#-i ../output/mm_cDNA_introns_97.idx \
#-o ../output/neuron10k_velocity \
#-x 10xv3 \
#-t 8 \
#neuron_10k_v3_S1_L002_R1_001.fastq.gz neuron_10k_v3_S1_L002_R2_001.fastq.gz \
#neuron_10k_v3_S1_L001_R1_001.fastq.gz neuron_10k_v3_S1_L001_R2_001.fastq.gz

date