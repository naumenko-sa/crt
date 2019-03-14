#!/bin/bash

ls *.junctions.txt | sed s/.junctions.txt// > table.txt

bam=`head -n1 table.txt`
prev_id=`qsub ~/crt/crt.load_junctions.localhd.pbs -v bam=$bam`

for bam in `cat table.txt | sed 1d`;
do
    job_id=`qsub ~/crt/crt.load_junctions.localhd.pbs -W depend=afterok:$prev_id -v bam=$bam`
    echo $bam $prev_id $job_id
    prev_id=$job_id
done
