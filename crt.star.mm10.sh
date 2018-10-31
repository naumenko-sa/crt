#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g

#modified from bcbio rna-seq log and STAR manual
#https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# len = mate lengths - 1

len=`zcat ${sample}_1.fq.gz | head -n 2| tail -n1 | awk '{print length-1}'`

echo "Overhang length:"$len

STAR --genomeDir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Mmusculus/mm10/star/ \
    --readFilesIn ${sample}_1.fq.gz ${sample}_2.fq.gz \
    --twopassMode Basic \
     --runThreadN 10 \
     --outFileNamePrefix $sample \
     --outReadsUnmapped Fastx \
     --outFilterMultimapNmax 10 \
     --outStd SAM  \
     --outSAMunmapped Within \
     --outSAMattributes NH HI NM MD AS  \
     --sjdbGTFfile /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf \
     --sjdbOverhang $len  --readFilesCommand zcat  --outSAMattrRGline ID:$sample PL:illumina PU:$sample SM:ID  \
     | samtools sort -@ 5 -m 1G  -T . \
     -o $sample.star.bam /dev/stdin

picard -Xmx20g MarkDuplicates \
    I=$sample.star.bam \
    O=$sample.bam \
    M=marked_dup_metrics.txt

samtools index $sample.bam

rm $sample.star.bam
