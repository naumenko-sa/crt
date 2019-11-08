#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=kenneth          # Job name
#SBATCH -c 5			    # cores
#SBATCH --mem-per-cpu=15G           # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date
echo $1

export LC_ALL=en_US.utf8 && export LANG=en_US.utf8 &&  /n/app/bcbio/dev/anaconda/bin/python /n/app/bcbio/dev/anaconda/bin/umis fasttagcount \
--cb_cutoff 1000 --genemap /n/scratch2/hsph_bioinformatic_core/sn240/no_trimming/ressler_mccullough2019/work/annotation/ref-transcripts-tx2gene.tsv  \
--cb_histogram /n/scratch2/hsph_bioinformatic_core/sn240/no_trimming/ressler_mccullough2019/work/umis/resslerkm-${1}/cb-histogram.txt \
--umi_matrix resslerkm-${1}-dupes.mtx.full  \
/n/scratch2/hsph_bioinformatic_core/sn240/no_trimming/ressler_mccullough2019/work/rapmap/resslerkm-${1}/resslerkm-${1}.bam \
resslerkm-${1}.mtx.full

date
export LC_ALL=en_US.utf8 && export LANG=en_US.utf8 &&  /n/app/bcbio/dev/anaconda/bin/python /n/app/bcbio/dev/anaconda/bin/umis sparse \
resslerkm-${1}.mtx.full resslerkm-${1}.mtx

export LC_ALL=en_US.utf8 && export LANG=en_US.utf8 &&  /n/app/bcbio/dev/anaconda/bin/python /n/app/bcbio/dev/anaconda/bin/umis sparse \
resslerkm-${1}-dupes.mtx.full resslerkm-${1}-dupes.mtx
