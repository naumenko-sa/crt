#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=5-00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=count            # Job name
#SBATCH -c 20			    # cores
#SBATCH --mem-per-cpu=3G            # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

date

ml cellranger/6.0.0

cellranger count \
--sample=$1 \
--id=$1 \
--fastqs=/path/to/03_mkfastq/ \
--transcriptome=/path/to/data/02_reference/refdata-gex-GRCh38-and-mm10-2020-A \
--localcores=20 \
--include-introns

date
