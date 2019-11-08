#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=kenneth       # Job name
#SBATCH -c 5			    # cores
#SBATCH --mem-per-cpu=15G           # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)


export LC_ALL=en_US.utf8 && export LANG=en_US.utf8 &&  /n/app/bcbio/dev/anaconda/bin/python \
/n/app/bcbio/tools/bin/umis demultiplex_samples --nedit 1 \
--barcodes /n/scratch2/hsph_bioinformatic_core/sn240/no_trimming/ressler_mccullough2019/work/umis/resslerkm/barcodes_missing.csv \
--out_dir /n/scratch2/hsph_bioinformatic_core/sn240/no_trimming/ressler_mccullough2019/work/umis/resslerkm/demultiplexed/barcodes_missing \
/n/scratch2/hsph_bioinformatic_core/sn240/no_trimming/ressler_mccullough2019/work/umis/resslerkm.umitransformed.fq.gz
