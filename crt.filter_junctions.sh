#!/bin/bash

module load python/3.5.2
module load sqlite/3.20.0

python3 ~/crt/FilterSpliceJunctions.py --sample $1 5 0.05

for f in *gtex_100;
do 
    cat $f | grep -v NULL > $f.no_NULL;
done;
