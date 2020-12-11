#!/bin/bash

# generate premrna gtf reference

bname=`basename $1 .gtf`

awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' $1 > $bname.premrna.gtf

# gffread -g /n/shared_db/bcbio/biodata/genomes/Sharrisii/msarhar1.11/seq/msarhar1.11.fa ref-transcripts.premrna.gtf -w msarhar1.11.unspliced.fa
