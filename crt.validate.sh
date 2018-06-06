#!/bin/bash

# validate vcf file against genome in a bottle calls for NA12878
# based on bcbio logs
# $1 = test.vcf.gz
# rtg manual
# https://github.com/RealTimeGenomics/rtg-tools/blob/master/installer/resources/tools/RTGOperationsManual.pdf

vpath=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878

bname=`basename $1 .vcf.gz`

gunzip -c $1 | grep "^#"  > $bname.nornaedit.vcf
gunzip -c $1 | grep -v "^#" | grep PASS | grep -v possible_rnaedit  >> $bname.nornaedit.vcf

bgzip $bname.nornaedit.vcf
tabix $bname.nornaedit.vcf

#uses PASS variants only
  export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=9100m && \
   rtg vcfeval --threads 5 \
   -b $vpath/truth_small_variants.vcf.gz \
   --bed-regions $vpath/truth_regions.bed \
   -c $bname.nornaedit.vcf.gz \
   -t /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
   -o rtg --vcf-score-field='GQ' 
#   --all-records

for f in {tp-baseline,fp,fn};
do
    echo snp $f `bcftools view --types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $bname.stat
    echo indels $f `bcftools view --exclude-types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $bname.stat
done


