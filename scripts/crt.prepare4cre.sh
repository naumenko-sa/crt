#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

# 50g is crucial -20,30 crashes sometimes

# prepares bcbio rna-seq variant for cre report generation
# case = family = project = S11 (example)
# cleanaup = 1

# qsub ~/crt/crt.prepare4cre.sh -v case=<project>,cleanup=1

function f_cleanup
{
    # better to look for project-summary than hardcode the year
    # keep bam files for new samples
    if [ -z $case ] 
    then
       echo "Project (case,family) folder does not exist. Exiting"
       exit 1
    fi

    cd $case
    
    project_summary=`find final -name project-summary.yaml`
    echo "Project summary: " $project_summary

    # if summary exists
    if [ -f $project_summary ] 
    then
        result_dir=`echo $project_summary | sed s/"\/project-summary.yaml"//`
        echo "Result dir:" $result_dir
    
        # if result_dir is empty that might cause copying entire /
        if [ -d $result_dir ] && [ -n "$result_dir"]
        then
            mv final/*/* .
	
	       cd ..
	       rm -rf ${case}/work
	       rm -rf ${case}/final
           rm -rf ${case}/transcriptome
	       cd $case

	       #rename bam files to match sample names
	       for f in *ready.bam;do mv $f `echo $f | sed s/"-ready"//`;done;
	       for f in *ready.bam.bai;do mv $f `echo $f | sed s/"-ready"//`;done;

	       rm ${sample}-transcriptome.bam
	            
	       #make bam files read only
	       for f in *.bam;do chmod 444 $f;done;
	        	        
	       #calculate md5 sums
	       for f in *.bam;do md5sum $f > $f.md5;done;

	       #validate bam files
	       for f in *.bam;do	cre.bam.validate.sh $f;done;
        fi
    fi
    
    cd ..
}

function f_get_sample_name
{
    cd $case
    sample=`cat samples.txt | head -n1`
    echo "Sample: " $sample
    export sample=$sample
    cd ..
}

function f_prepare
{
    #downstream scripts use $vcf for current vcf file
    cd $case
    original_vcf=${case}-gatk-haplotype-annotated.vcf.gz

    # to run cre immediately after that
    mkdir $case
    
    cp $original_vcf ${original_vcf}.tbi $case

    cd $case
    
    gunzip -c $original_vcf | grep "^#"  > $sample.vcf
    gunzip -c $original_vcf | grep -v "^#" | grep PASS | grep -v possible_rnaedit | egrep -v "^GL000" >> $sample.vcf

    bgzip $sample.vcf
    tabix $sample.vcf.gz

    cre.vt.decompose.sh $sample.vcf.gz
    cre.vep.sh $sample.decomposed.vcf.gz
    cre.vcfanno.sh $sample.decomposed.vepeffects.vcf.gz
    echo -e $case"\t"$sample"\t0\t0\t0\t0\n" > $case.ped
    vcf2db.py $sample.decomposed.vepeffects.annotated.vcf.gz 	$case.ped	${case}-ensemble.db

    mv $sample.decomposed.vepeffects.annotated.vcf.gz ${case}-ensemble-annotated-decomposed.vcf.gz
    mv $sample.decomposed.vepeffects.annotated.vcf.gz.tbi ${case}-ensemble-annotated-decomposed.vcf.gz.tbi

    ln -s ${case}-ensemble-annotated-decomposed.vcf.gz ${case}-gatk-haplotype-annotated-decomposed.vcf.gz
    ln -s ${case}-ensemble-annotated-decomposed.vcf.gz.tbi ${case}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi
    
    cd ..
    cd ..
}

if [ -z $case ] || [ ! -d $case ]
then
    echo "Case (project,family) folder does not exist. Exiting"
    exit 1
else
    f_get_sample_name
    if [ $cleanup -eq 1 ]
    then
        f_cleanup
    fi
    
    f_prepare
fi
