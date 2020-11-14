#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=short           # Partition (queue) priority
#SBATCH --time=10:00:00             # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=barcode_stats    # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=10G                   # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

# calculate barcode csv table

# $1 = fq.gz

gunzip -c $1 | awk '{if (NR%4==2) print $0}' | sort | uniq -c | awk '{print $2","$1}' | sort -t "," -k2,2nr > $1.stats.csv
