#!/bin/bash

PHENO=$1
ANCESTRY=$2
CHROM=$3
CHUNKS_DIR=$4
OUT_DIR="."

if [[ $CHROM == 23 ]]; then
    CHROM="X"
fi

PATTERN="*chr${CHROM}_*_${ANCESTRY}.psam"
find $CHUNKS_DIR -name $PATTERN | sed 's/\.[^.]*$//' > mergelist_${CHROM}_${ANCESTRY}.txt

/opt/conda/bin/plink2 \
    --pmerge-list mergelist_${chrom}_${ANCESTRY}.txt \
    --out ${OUT_DIR}/${PHENO}_chr${chrom}_${ANCESTRY}