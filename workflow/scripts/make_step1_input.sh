#!/bin/bash

PHENO=$1
ANCESTRY=$2
CHROM_DIR=$3
OUT_DIR="."

chromlist=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
for chrom in "${chromlist[@]}"; do
    /opt/conda/bin/plink2 \
        --pfile ${CHROM_DIR}/${PHENO}_chr${chrom}_${ANCESTRY} \
        --indep-pairwise 50 5 0.05 \
        --out ${OUT_DIR}/${PHENO}_chr${chrom}_${ANCESTRY} \

    /opt/conda/bin/plink2 \
        --pfile ${CHROM_DIR}/${PHENO}_chr${chrom}_${ANCESTRY} \
        --extract ${OUT_DIR}/${PHENO}_chr${chrom}_${ANCESTRY}.prune.in \
        --make-pgen \
        --out ${OUT_DIR}/${PHENO}_chr${chrom}_${ANCESTRY}   
done

find ${OUT_DIR}/ -maxdepth 1 -name "${PHENO}_*_${ANCESTRY}.psam" | sed 's/\.[^.]*$//' > pruned_files_${ANCESTRY}.txt
/opt/conda/bin/plink2 \
    --pmerge-list pruned_files_${ANCESTRY}.txt \
    --make-bed \
    --out ${OUT_DIR}/${PHENO}_${ANCESTRY}_pruned
