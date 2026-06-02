#!/bin/bash

PHENO=$1
ANCESTRY=$2
PHENO_PATH=$3
COVARS=$4
QCOVARS=$5
PLINK_PATH=$6
TRAIT_TYPE=$7
INV_NORMALIZE=$8
N_CPUS=$9
OUT_DIR="."

/usr/local/bin/step1_fitNULLGLMM.R \
	--plinkFile=${PLINK_PATH} \
    --useSparseGRMtoFitNULL=FALSE \
    --phenoFile=${PHENO_PATH} \
    --phenoCol=${PHENO} \
    --covarColList=${COVARS} \
    --qCovarColList=${QCOVARS}  \
    --sampleIDColinphenoFile=IID \
    --invNormalize=FALSE \
    --traitType=quantitative        \
    --outputPrefix=$OUT_DIR/saige_step1_${PHENO}_${ANCESTRY} \
    --nThreads=$N_CPUS \
    --IsOverwriteVarianceRatioFile=TRUE