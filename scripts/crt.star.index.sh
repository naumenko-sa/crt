#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=5-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=schumacher            # Job name
#SBATCH -c 20
#SBATCH --mem=30G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

# mouse input files: https://www.gencodegenes.org/mouse/

STAR \
--genomeDir /n/data1/cores/bcbio/PIs/valerie_schumacher/schumacher_RNAediting_analysis_in_mouse_podocytes_after_kidney_injury_hbc03712/2019-08-08_star_from_scratch \
--genomeFastaFiles /n/data1/cores/bcbio/PIs/valerie_schumacher/schumacher_RNAediting_analysis_in_mouse_podocytes_after_kidney_injury_hbc03712/2019-08-08_star_from_scratch/GRCm38.primary_assembly.genome.fa \
--runThreadN 20 \
--limitGenomeGenerateRAM 30000000000 \
--genomeChrBinNbits 14 \
--runMode genomeGenerate \
--genomeSAindexNbases 14 \
--sjdbOverhang 149 \
--sjdbGTFfile /n/data1/cores/bcbio/PIs/valerie_schumacher/schumacher_RNAediting_analysis_in_mouse_podocytes_after_kidney_injury_hbc03712/2019-08-08_star_from_scratch/gencode.vM22.annotation.gtf
