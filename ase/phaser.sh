#!/bin/bash

# https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/

date

module load gcc/6.2.0
module load python/2.7.12
module load cython/0.25.1

which python

export PHASER_PATH=/n/data1/cores/bcbio/naumenko/tools/phaser

python $PHASER_PATH/phaser/phaser.py \
--vcf NA06986.vcf.gz \
--bam NA06986.2.M_111215_4.bam \
--paired_end 1 \
--mapq 255 \
--baseq 10 \
--sample NA06986 \
--blacklist hg19_hla.bed \
--haplo_count_blacklist hg19_haplo_count_blacklist.bed \
--threads 1 --o phaser_test_case

python $PHASER_PATH/phaser_gene_ae/phaser_gene_ae.py \
--haplotypic_counts phaser_test_case.haplotypic_counts.txt \
--features gencode.v19.GRCh37.genes.bed \
--o phaser_test_case_gene_ae.txt