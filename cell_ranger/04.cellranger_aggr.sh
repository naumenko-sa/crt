#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=5-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=aggr             # Job name
#SBATCH -c 20			    # cores
#SBATCH --mem-per-cpu=3G            # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

date

ml cellranger/6.0.0

cellranger aggr \
--csv=/n/data1/cores/bcbio/PIs/rachel_clark/hbc_10x_scRNAseq_Clark_gy_irradition_skin_xenografts_human_hbc04213/data/05_cellranger_aggr/samples.csv \
--id=M_GY_0_2 \
--localcores=10

date
