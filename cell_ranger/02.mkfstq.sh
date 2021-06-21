#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=10:00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=mkfastq          # Job name
#SBATCH -c 10			    # cores
#SBATCH --mem-per-cpu=3G            # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

date

ml cellranger/6.0.0
ml bcl2fastq/2.20.0.422

cellranger mkfastq \
--id=clark2021 \
--run=/path/to/01_raw_data \
--csv=/path/to/meta/meta.csv \
--localcores=10 \
--output-dir=/path/to/data/03_mkfastq

date
