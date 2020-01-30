#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=medium        # Partition (queue)
#SBATCH --time=5-00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=sc_step2         # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=50G            # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

# demultiplex samples
# 4.5bln reads = 42h crashed with 50G ram, passed w100G ram
# 1.4 bln = 10h?
# 2.2 bln = 18h
# 1 bln = 8h20 min

date
export file_fq=$1
export barcodes=$2
export LC_ALL=en_US.utf8 && export LANG=en_US.utf8 && \
/n/app/bcbio/dev/anaconda/bin/python /n/app/bcbio/dev/anaconda/bin/umis demultiplex_samples --nedit 1 \
--barcodes $barcodes \
--out_dir demultiplexed \
$file_fq

date

