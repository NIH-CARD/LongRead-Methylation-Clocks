#!/bin/bash

COHORT=$1
WORK_DIR=$2
RESULTS_DIR="${WORK_DIR}/results"
DATA_DIR="${RESULTS_DIR}/preprocessing"


suffixes=("promoter" "promoter_under80" "promoter_clock1" "promoter_clock2" "promoter_clock3" "CpG_clock1" "CpG_clock2" "CpG_clock3")

for suffix in "${suffixes[@]}"; do
    genoml continuous supervised munge --prefix ${RESULTS_DIR}/${COHORT}_${suffix} --addit ${DATA_DIR}/${COHORT}/${COHORT}_train_${suffix}.tsv --pheno ${DATA_DIR}/${COHORT}/${COHORT}_pheno.tsv --addit_test ${DATA_DIR}/${COHORT}/${COHORT}_test_${suffix}.tsv --pheno_test ${DATA_DIR}/${COHORT}/${COHORT}_pheno.tsv &&
    genoml continuous supervised train --prefix ${RESULTS_DIR}/${COHORT}_${suffix} --verbose &&
    genoml continuous supervised tune --prefix ${RESULTS_DIR}/${COHORT}_${suffix} --verbose &&
    genoml continuous supervised test --prefix ${RESULTS_DIR}/${COHORT}_${suffix} --verbose
done
