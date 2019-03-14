#!/bin/bash
# extract variants in protein coding genes for allelic imbalance analysis

if [ -z $vcf ]
then
    vcf=$1
fi

bname=`basename $vcf .vcf.gz`

echo "chrom,pos,ref,alt" > $bname.het.csv

bcftools view -v snps -R ~/cre/data/protein_coding_genes.exons.bed $vcf | \
bcftools query -i 'FMT/GT="0/1" && INFO/DP>=20 && INFO/possible_rnaedit=0 && FILTER="PASS"' -f '%CHROM,%POS,[%AD]\n' - >> $bname.het.csv
