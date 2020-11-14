#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=umitransform     # Job name
#SBATCH -c 7			    # cores
#SBATCH --mem=50G                   # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date

# this step is needed before separation of multiseq and transcript libraries

umis_prefix=/n/app/bcbio/dev/anaconda/lib/python3.6/site-packages/bcbio/data/umis

umis fastqtransform \
--separate_cb $umis_prefix/harvard-indrop-v3-transform.json $1 $2 $3 $4 | gzip > umitransformed.fq.gz
