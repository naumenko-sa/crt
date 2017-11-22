#!/bin/bash
#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

#50g is crucial -20,30 crashes sometimes

#prepare bcbio rna-seq variant to cre report generation
#put 

if [ -z $sample ]
then
    sample=$1
fi

#downstream scripts use $vcf for current vcf file
if [ -z $original_vcf ]
then
    original_vcf=$2
fi

if [ -z $family ]
then
    family=$3
fi

gunzip -c $original_vcf | grep "^#"  > $sample.vcf
gunzip -c $original_vcf | grep -v "^#" | grep PASS | grep -v possible_rnaedit  >> $sample.vcf

bgzip $sample.vcf
tabix $sample.vcf.gz

cre.vt.decompose.sh $sample.vcf.gz
cre.vep.sh $sample.decomposed.vcf.gz
cre.gemini_load.sh $sample.decomposed.vepeffects.vcf.gz
#gemini.gemini2txt.sh $sample.decomposed.vepeffects.db

mv $sample.decomposed.vepeffects.db ${family}-ensemble.db
mv $sample.decomposed.vepeffects.vcf.gz ${family}-ensemble-annotated-decomposed.vcf.gz
mv $sample.decomposed.vepeffects.vcf.gz.tbi ${family}-ensemble-annotated-decomposed.vcf.gz.tbi

ln -s ${family}-ensemble-annotated-decomposed.vcf.gz ${family}-gatk-haplotype-annotated-decomposed.vcf.gz
ln -s ${family}-ensemble-annotated-decomposed.vcf.gz.tbi ${family}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi
