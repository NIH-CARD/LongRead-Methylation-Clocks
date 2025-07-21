#!/bin/bash

COHORT=$1
WORK_DIR=$2
DATA_DIR=$3
RESULTS_DIR=${WORK_DIR}/results
PREPROCESSING_DIR=${WORK_DIR}/results/preprocessing

mkdir ${RESULTS_DIR}/${COHORT}_promoter
mkdir ${RESULTS_DIR}/${COHORT}_promoter_clock1
mkdir ${RESULTS_DIR}/${COHORT}_promoter_clock2
mkdir ${RESULTS_DIR}/${COHORT}_promoter_clock3
mkdir ${RESULTS_DIR}/${COHORT}_CpG_clock1
mkdir ${RESULTS_DIR}/${COHORT}_CpG_clock2
mkdir ${RESULTS_DIR}/${COHORT}_CpG_clock3
mkdir ${RESULTS_DIR}/${COHORT}_promoter_under80

mkdir ${PREPROCESSING_DIR}/${COHORT}
mkdir ${PREPROCESSING_DIR}/${COHORT}/promoter_data
mkdir ${PREPROCESSING_DIR}/${COHORT}/promoter_data/chunks
mkdir ${PREPROCESSING_DIR}/${COHORT}/promoter_data/sig_chunks
mkdir ${PREPROCESSING_DIR}/${COHORT}/promoter_data_under80
mkdir ${PREPROCESSING_DIR}/${COHORT}/promoter_data_under80/chunks
mkdir ${PREPROCESSING_DIR}/${COHORT}/promoter_data_under80/sig_chunks


# HBCC data are phased and include duplicate CpGs, while NABEC are unphased and do not contain duplicates
if [[ "${COHORT}" == "HBCC" ]]; then
    for hap in {1..2}; do
        gunzip -c ${DATA_DIR}/HBCC_combined_methylation_hap${hap}.042025.tsv.sorted.tsv.gz.Hs_EPDnew_006_hg38.2kb.6cols.bed.tsv.gz \
        | awk '!seen[$0]++' \
        > ${PREPROCESSING_DIR}/HBCC/HBCC_CpGs_hap${hap}.tsv 
    done
elif [[ "${COHORT}" == "NABEC" ]]; then
    gunzip -c ${DATA_DIR}/NABEC_combinedHaps_06182025.Hs_EPDnew_006_hg38.2kb.6cols.sorted.tsv.gz > ${PREPROCESSING_DIR}/NABEC/NABEC_CpGs_merged.tsv
fi
