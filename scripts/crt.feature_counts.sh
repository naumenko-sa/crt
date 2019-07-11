#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# calculates features (reads) for RPKM calculation in R, outputs length of the genes
# we were using RPKMs not FPKMs to be in line with Cummings2017 and Kremer2017.
# in bcbio featureCounts counts fragments

# bam=file.bam
# strandedness = [0,1]

# by default it counts exon features, use 
# -t transcript or -t three_prime_utr
# for custom gtf

if [ -z $bam ]
then
    bam=$1
fi

if [ -z $gtf ]
then
    gtf=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf
fi

if [ -z $strandedness ]
then
    strandedness=0
fi

echo "bam: " $bam
echo "strandedness: " $strandedness

# -s 0 unstranded (
# -p - count pairs (fragments) not reads
# -B only properly paired (both ends mapped)
# -C don't count mapped to different chromosomes or different strands
featureCounts -T 8 \
	-t three_prime_utr \
	-g gene_id \
	-s $strandedness \
	-C \
	--largestOverlap \
	-a $gtf \
	-o $bam.feature_counts.txt $bam
	#-s 0 -p -B -C
