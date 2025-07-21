#!/bin/bash

WORK_DIR=${1}
DATA_DIR=${2}
COV_THRESHOLD_READS=${3}
COV_THRESHOLD_RATIO=${4}
N_PCS=${5}
CORR_CHUNK_SIZE=${6}
PEARSON_CORR_THRESHOLD=${7}
TARGET_COLUMN=${8}
RUN_VIF=${9}
VIF_THRESHOLD=${10}
PROMOTER_LENGTH=${11}
NUM_SHAP_VALS=${12}

COHORTS=("NABEC" "HBCC")
ANALYSIS_DIR="${WORK_DIR}/analyses"

${ANALYSIS_DIR}/00_initialize_directories.sh
for COHORT in "${COHORTS[@]}"; do
    ${ANALYSIS_DIR}/01_unzip_files.sh $COHORT $WORK_DIR $DATA_DIR
    ${ANALYSIS_DIR}/02_merge_haplotypes.py --cohort $COHORT --work_dir $WORK_DIR
    ${ANALYSIS_DIR}/03_aggregate_promoters.py --cohort $COHORT --work_dir $WORK_DIR --data_dir $DATA_DIR --cov_reads $COV_THRESHOLD_READS --cov_ratio $COV_THRESHOLD_RATIO
    ${ANALYSIS_DIR}/04_calculate_pcs.py --cohort $COHORT --work_dir $WORK_DIR --n_pcs $N_PCS
    ${ANALYSIS_DIR}/05_process_metadata.py --cohort $COHORT --work_dir $WORK_DIR
    ${ANALYSIS_DIR}/06_metadata_promoter_merge.py --cohort $COHORT --work_dir $WORK_DIR
    ${ANALYSIS_DIR}/07_promoter_regressions.py --cohort $COHORT --work_dir $WORK_DIR --chunk_size $CORR_CHUNK_SIZE --corr_threshold $PEARSON_CORR_THRESHOLD --target_column $TARGET_COLUMN --vif $RUN_VIF --vif_threshold $VIF_THRESHOLD
    ${ANALYSIS_DIR}/08_merge_sig_promoters.py --cohort $COHORT --work_dir $WORK_DIR
    ${ANALYSIS_DIR}/09_extract_clocks.py --cohort $COHORT --work_dir $WORK_DIR --data_dir $DATA_DIR --promoter_length $PROMOTER_LENGTH
    ${ANALYSIS_DIR}/10_prep_for_genoml.py --cohort $COHORT --work_dir $WORK_DIR
    ${ANALYSIS_DIR}/11_exec_genoml.sh $COHORT $WORK_DIR
    ${ANALYSIS_DIR}/12_significant_genes.py --cohort $COHORT --work_dir $WORK_DIR
    ${ANALYSIS_DIR}/13_shap.py --cohort $COHORT --work_dir $WORK_DIR --n_shap $NUM_SHAP_VALS
done
${ANALYSIS_DIR}/14_gene_promoter_overlaps.py --work_dir $WORK_DIR
${ANALYSIS_DIR}/15_generate_figs_and_tables.py --work_dir $WORK_DIR
