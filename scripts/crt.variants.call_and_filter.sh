#!/bin/bash
# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=short           # Partition (queue) priority
#SBATCH --time=10:00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=gatk          # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=20G           # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

# call varians from split'n'cigar prepared bam
# input:
# $1 - input.bam
# output:
# input.vcf.gz

bname=`basename $1 .bam`

reference=/n/shared_db/bcbio/biodata/genomes/Hsapiens/hg38/seq/hg38.fa

gatk3 -Xmx20g -T HaplotypeCaller \
   -R $reference \
   -I $1 \
   -dontUseSoftClippedBases \
   -stand_call_conf 20 \
   -o $bname.raw.vcf \
   -L chrX \
   -L chrY

gatk3 -Xmx20g -T VariantFiltration \
    -R $reference \
    -V $bname.raw.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" \
    -filterName QD -filter "QD < 2.0" \
    -o $bname.vcf

bgzip $bname.vcf
tabix $bname.vcf.gz
