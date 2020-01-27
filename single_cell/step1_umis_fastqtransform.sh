#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=medium        # Partition (queue)
#SBATCH --time=5-00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=sc_step1         # Job name
#SBATCH -c 16			    # cores
#SBATCH --mem-per-cpu=3G            # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

# shift barcodes to read name
# 4.5bln reads = 41h
# 2.2bln reads = 30h
# 9991mln reads = 12h

date
export sample=$1
export LC_ALL=en_US.utf8 && export LANG=en_US.utf8 &&  \
/n/app/bcbio/dev/anaconda/bin/python \
/n/app/bcbio/tools/bin/umis fastqtransform \
--separate_cb /n/app/bcbio/dev/anaconda/lib/python3.6/site-packages/bcbio_nextgen-1.2.0a0-py3.6.egg/bcbio/rnaseq/../data/umis/harvard-indrop-v3-transform.json \
--cores 16  \
${sample}_1.fq.gz \
${sample}_2.fq.gz \
${sample}_3.fq.gz \
${sample}_4.fq.gz | seqtk seq -L 20 - | gzip > ${sample}.umitransformed.fq.gz

date

