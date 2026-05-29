#!/bin/bash

chromlist=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 "X")

FILE_DIR=$1
ANCESTRY=$2

mkdir -p mergelists

for chrom in "${chromlist[@]}"
do
    echo "Writing chr${chrom} to mergelist_${chrom}_${ANCESTRY}.txt"
    PATTERN="*chr${chrom}_*_${ANCESTRY}.psam"
    find $FILE_DIR -name $PATTERN | sed 's/\.[^.]*$//' > mergelists/mergelist_${chrom}_${ANCESTRY}.txt
done
