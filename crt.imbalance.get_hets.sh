#!/bin/bash
# extracts variants in protein coding genes for allelic imbalance analysis

echo "chrom,pos,ref,alt" > $1.het.csv

bcftools view -v snps -R ~/cre/data/protein_coding_genes.exons.bed $1 | \
bcftools query -i 'FMT/GT="0/1" && INFO/DP>=20 && INFO/possible_rnaedit=0 && FILTER="PASS"' -f '%CHROM,%POS,[%AD]\n' - >> $1.het.csv
