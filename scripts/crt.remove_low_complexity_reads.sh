#!/bin/bash

#PBS -l walltime=15:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=30g,mem=30g

# params:
# left, right, output

~/crt/prinseq-lite.pl -fastq $left -fastq2 $right -log $output.prinseq.log -derep 12345 -lc_method dust -lc_threshold 7 -out_good $output