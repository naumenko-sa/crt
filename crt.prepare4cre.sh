#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

#50g is crucial -20,30 crashes sometimes

# prepares bcbio rna-seq variant for cre report generation
# $1 = family = project = case = S11 (example)


function f_cleanup
{
    # better to look for project-summary than hardcode the year
    # keep bam files for new samples
    cd $family
    result_dir=`find final -name project-summary.yaml | sed s/"\/project-summary.yaml"//`
    echo $result_dir
    
    #just to be on the safe side
    if [ -d $result_dir ]
    then
	mv final/*/* .
	rmdir final
	
	cd ..
	rm -rf ${family}/work
	cd $family
	
	#rename bam files to match sample names
	for f in *ready.bam;do mv $f `echo $f | sed s/"-ready"//`;done;
	for f in *ready.bam.bai;do mv $f `echo $f | sed s/"-ready"//`;done;
	            
	#make bam files read only
	for f in *.bam;do chmod 444 $f;done;
	        	        
	#calculate md5 sums
	for f in *.bam;do md5sum $f > $f.md5;done;

	#validate bam files
	for f in *.bam;do	cre.bam.validate.sh $f;done;
    fi
    
    cd ..
}


function f_prepare
{

    #downstream scripts use $vcf for current vcf file
    cd $family
    sample=`cat samples.txt | head -n1`
    echo "Sample: " $sample
    original_vcf=${family}-gatk-haplotype-annotated.vcf.gz

    # to run cre immediately after that
    mkdir $family
    
    cp $original_vcf ${original_vcf}.tbi $family

    cd $family
    
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
    
    cd ..
    cd ..
}

if [ -z $family ] || [ ! -d $family ]
then
    echo "Project (family) folder does not exist. Exiting"
    exit 1
else
    if [ $cleanup -eq 1 ]
    then
        f_cleanup
    fi
    
    f_prepare
fi