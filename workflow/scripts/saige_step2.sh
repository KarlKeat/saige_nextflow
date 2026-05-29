#!/bin/bash

PHENO=$1
ANCESTRY=$2
PHENO_PATH=$3
GRM_FILE=$4
VARIANCE_RATIO_FILE=$5
PGEN_DIR=$6
OUT_DIR="."

chromlist=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
for chrom in "${chromlist[@]}"; do
    /usr/local/bin/step2_SPAtests.R \
        --pgenPrefix=${PGEN_DIR}/${PHENO}_pmbb_chr${chrom}_${ANCESTRY} \
        --sampleFile=${PHENO_PATH} \
        --AlleleOrder=ref-first \
        --SAIGEOutputFile=${OUT_DIR}/saige_step2_${PHENO}_${ANCESTRY}_chr${chrom}.txt \
        --chrom=$chrom \
        --minMAF=0 \
        --minMAC=20 \
        --GMMATmodelFile=${GRM_FILE} \
        --varianceRatioFile=${VARIANCE_RATIO_FILE} \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --is_output_moreDetails=TRUE \
        --LOCO=TRUE
done