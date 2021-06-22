#!/bin/bash

cat tagcounts.mtx | sed 1,3d | awk '{sum+=$3}END{print sum}' > tagcounts.mtx.total_counts
cat tagcounts-dupes.mtx | sed 1,3d | awk '{sum+=$3}END{print sum}' > tagcounts-dupes.mtx.total_counts

