#!/bin/bash

PLINK_PREFIX=$1
ANCESTRY=$2
PHENO_COVAR_FILE=$3
OUT_DIR="."

BASENAME="${PLINK_PREFIX##*/}"

cut -f2 ${PHENO_COVAR_FILE} > temp_ids.txt

/opt/conda/bin/plink2 \
	--pfile ${PLINK_PREFIX} \
	--maf 0.01 \
	--geno 0.05 \
	--keep temp_ids.txt \
	--hwe 1e-6 keep-fewhet 0.0001 \
	--make-pgen \
	--out ${OUT_DIR}/${BASENAME}_filtered_${ANCESTRY}

