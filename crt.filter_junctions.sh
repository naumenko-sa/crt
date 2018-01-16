#!/bin/bash

module load python/3.5.2
module load sqlite/3.20.0

python3 ~/crt/FilterSpliceJunctions.py --sample $1 5 0.05

