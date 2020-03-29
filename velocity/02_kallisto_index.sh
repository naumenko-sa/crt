#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=kallisto          # Job name
#SBATCH -c 8			    # cores
#SBATCH --mem=100G           # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date
which kallisto

# killed with 60G of ram
# passed with 100G

# https://bustools.github.io/BUS_notebooks_R/velocity.html

kallisto index -i mm_cDNA_introns_97.idx cDNA_introns.fa

date
