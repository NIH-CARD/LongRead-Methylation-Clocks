#!/bin/bash

WORK_DIR="/data/CARDPB/projects/METH_MODEL"
RESULTS_DIR="${WORK_DIR}/results"
ANALYSIS_DIR="${WORK_DIR}/analyses"
PREPROCESSING_DIR="${WORK_DIR}/results/preprocessing"

mkdir -p ${PREPROCESSING_DIR}

wget -P ${PREPROCESSING_DIR} https://epd.expasy.org/ftp/epdnew/H_sapiens/current/Hs_EPDnew.dat
grep "^ID" ${PREPROCESSING_DIR}/Hs_EPDnew.dat | awk '{print $2}' > ${PREPROCESSING_DIR}/promoter_ids.txt
grep "^GN" ${PREPROCESSING_DIR}/Hs_EPDnew.dat | awk -F'[=;]' '{print $2}' > ${PREPROCESSING_DIR}/promoter_gene_names.txt
paste ${PREPROCESSING_DIR}/promoter_ids.txt ${PREPROCESSING_DIR}/promoter_gene_names.txt > ${PREPROCESSING_DIR}/promoter_gene_mapping.txt

rm ${PREPROCESSING_DIR}/Hs_EPDnew.dat
rm ${PREPROCESSING_DIR}/promoter_ids.txt
rm ${PREPROCESSING_DIR}/promoter_gene_names.txt
