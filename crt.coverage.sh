#!/bin/bash
# get coverage for exons of a gene
# $sample - sample name, $sample.chrom.bam should be in the current directory
# $gene - gene, gene.bed should be in the current directory
# output - $sample.$gene.coverage
# first generate bam files for each chromosome:
# samtools view -bh DMD.bam X > DMD.X.bam
# samtools index DMD.X.bam
# otherwise (1chr for bed, 1chr for bam it is not working)
# -split for RNA-seq coverage

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

if [ -z $sample ]
then
    sample=$1
fi

if [ -z $bed ]
then
    bed=$2
fi

if [ -z $rnaseq ]
then
    rnaseq=
else
    rnaseq='-split'
fi

echo $rnaseq

chrom=`head -n1 $bed | cut -f1`

echo "Calculates coverage of" $bed " in " $sample " for chr " $chrom

echo -e "chrom\tstart\tend\texon\t$sample" > $sample.$bed.coverage
bedtools coverage -a $bed -b $sample.$chrom.bam $rnaseq -mean -sorted >> $sample.$bed.coverage
