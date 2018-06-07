#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# validate vcf file against genome in a bottle calls for NA12878
# based on bcbio logs
# $1 = test.vcf.gz
# rtg manual
# https://github.com/RealTimeGenomics/rtg-tools/blob/master/installer/resources/tools/RTGOperationsManual.pdf

if [ -z $vcf ]
then
    vcf=$1
fi

vpath=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878
protein_coding_exons=/home/naumenko/Desktop/reference_tables/protein_coding_genes.exons.fixed.bed

bname=`basename $1 .vcf.gz`

bcftools filter -e "INFO/DP<10 || possible_rnaedit=1 || FILTER!='PASS'" $1 > $bname.nornaediting.vcf
bgzip $bname.nornaediting.vcf
tabix $bname.nornaediting.vcf.gz

bedtools intersect -a $bname.nornaediting.vcf.gz -b $vpath/truth_regions.bed -header > $bname.in_truth_regions.vcf
bgzip $bname.in_truth_regions.vcf
tabix $bname.in_truth_regions.vcf.gz

bedtools intersect -a $bname.in_truth_regions.vcf.gz -b $protein_coding_exons -header > $bname.filtered.vcf
bgzip $bname.filtered.vcf
tabix $bname.filtered.vcf.gz

#uses PASS variants only
export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=9100m && \
   rtg vcfeval --threads 5 \
   -b $vpath/truth_small_variants.vcf.gz \
   --bed-regions $vpath/truth_regions.bed \
   -c $bname.filtered.vcf.gz \
   -t /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
   -o rtg --vcf-score-field='GQ' 
#   --all-records

for f in {tp-baseline,fp,fn};
do
    echo snp $f `bcftools view --types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $bname.stat
    echo indels $f `bcftools view --exclude-types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $bname.stat
done
