#!/bin/bash

#compares variants detected in trios

for f in *.panel.csv;
do
    cat $f | sed 1d | awk -F '","' '{if ($13>=10) print $1"-"$3"-"$4}' | sed s/\"//g | sort > $f.vars
done
