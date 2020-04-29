#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=medium        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=kallisto          # Job name
#SBATCH -c 8			    # cores
#SBATCH --mem=100G           # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date

source activate r
which kb

# https://bustools.github.io/BUS_notebooks_R/velocity.html

kreference_prefix=/n/data1/cores/bcbio/naumenko/velocity_test/veloindex_indrops3

# 100G ram
kb count \
-i ${kreference_prefix}/mm_cDNA_introns_97.idx \
-g ${kreference_prefix}/tr2g.tsv \
-x INDROPSV3 \
-o kallisto_bus_output \
-c1 ${kreference_prefix}/cDNA_tx_to_capture.txt \
-c2 ${kreference_prefix}/introns_tx_to_capture.txt \
--lamanno \
--verbose \
-t 8 \
${1}_2.fq.gz ${1}_4.fq.gz ${1}_1.fq.gz

source deactivate

date