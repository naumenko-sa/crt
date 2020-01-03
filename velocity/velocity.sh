#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=kallisto          # Job name
#SBATCH -c 8			    # cores
#SBATCH --mem=10G           # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date
which kallisto

# killed with 60G of ram
# passed with 100G

# https://bustools.github.io/BUS_notebooks_R/velocity.html

# kallisto index -i ./output/mm_cDNA_introns_97.idx ./output/neuron10k_velocity/cDNA_introns.fa
# 100G ram
#kallisto bus \
#-i ../output/mm_cDNA_introns_97.idx \
#-o ../output/neuron10k_velocity \
#-x 10xv3 \
#-t 8 \
#neuron_10k_v3_S1_L002_R1_001.fastq.gz neuron_10k_v3_S1_L002_R2_001.fastq.gz \
#neuron_10k_v3_S1_L001_R1_001.fastq.gz neuron_10k_v3_S1_L001_R2_001.fastq.gz

# step 
# barcode correction and capturing reads on the introns and on the exons
# 10 min
# bustools correct \
# -w whitelist_v3.txt \
# -p output.bus | \
# bustools sort -o output.correct.sort.bus -t8 -
# bustools capture -s -x -o spliced.bus -c introns_tx_to_capture.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
# bustools capture -s -x -o unspliced.bus -c cDNA_tx_to_capture.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus

# 3 min 
# 10G ram
bustools count -o unspliced -g tr2g.tsv -e matrix.ec -t transcripts.txt --genecounts unspliced.bus
bustools count -o spliced -g tr2g.tsv -e matrix.ec -t transcripts.txt --genecounts spliced.bus

date
