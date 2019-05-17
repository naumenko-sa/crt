#!/bin/bash

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

module load python/3.5.2
module load sqlite/3.20.0

if [ -z $bam ]
then
    bam=$1
fi

python3 ~/crt/scripts/FilterSpliceJunctions.py --sample $bam 5 0.05

cat ${bam}_specific_rc5_norm_rc0.05_n_gtex_184 | grep -v BOTH > $bam.rare_junctions.txt
