#!/bin/bash

#junction statistics

echo "Sample" `echo $1 | awk -F "." '{print $1}'`
echo "Start" `cat $1 | grep -c START`
echo "End" `cat $1 | grep -c STOP`
echo "Exon_skipping" `cat $1 | grep -c EXON_SKIP`
echo "BOTH" `cat $1 | grep -c BOTH`
echo "Total" `cat $1 | sed 1d | wc -l`
echo "Genes" `cat $1 | sed 1d | awk '{print $1}' | sort | uniq | wc -l`