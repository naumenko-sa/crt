#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=10:00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=kenneth          # Job name
#SBATCH -c 16			    # cores
#SBATCH --mem-per-cpu=3G            # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date
rapmap quasimap -t 16 \
-i /path/project/work/rapmap/index/quasiindex/mm10 \
-r /path/project/work/umis/resslerkm-${1}.filtered.fq.gz | \
samtools sort -@ 16 -m 1G  -o resslerkm-${1}.bam /dev/stdin
samtools index -@ 16 resslerkm-${1}.bam resslerkm-${1}.bam.bai
date
