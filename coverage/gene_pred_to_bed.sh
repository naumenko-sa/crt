#!/bin/bash

# convert gene pred format to bed usable by seqc

cat $1 | awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' > $1.bed
