#!/bin/bash

#fixes issue https://github.com/chapmanb/bcbio-nextgen/issues/2135

#bcbio rnaseq crashes

cd work/variation

mv combined.vcf.gz combined.before_fixing.vcf.gz
mv combined.vcf.gz.tbi combined.before_fixing.vcf.gz.tbi

bcftools view -o combined.vcf.gz -O z combined.before_fixing.vcf.gz 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y
tabix combined.vcf.gz

cd ../../

#rerun bcbio