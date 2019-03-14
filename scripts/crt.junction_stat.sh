#!/bin/bash

#junction statistics

cat $1 | grep -v NULL > $1.not_NULL

Sample=`echo $1.not_NULL | awk -F "." '{print $1}'`
#new start is where only stop (end) is annotated
Start=`cat $1.not_NULL | grep -c STOP`
#new junction end is where only start is annotated
End=`cat $1.not_NULL | grep -c START`
Exon_skipping=`cat $1.not_NULL | grep -c EXON_SKIP`
BOTH=`cat $1.not_NULL | grep -c BOTH`
Total=`cat $1.not_NULL | sed 1d | wc -l`
Genes=`cat $1.not_NULL | sed 1d | awk '{print $1}' | sort | uniq | wc -l`

echo -e "Sample\tStart\tEnd\tExon_skipping\tBOTH\tTotal\tGenes"
echo -e $Sample"\t"$Start"\t"$End"\t"$Exon_skipping"\t"$BOTH"\t"$Total"\t"$Genes
