#!/usr/bin/env python3


# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import numpy as np
import os
import pandas as pd
import re
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor


# =========================================================================
# PARSE PARAMS PASSED BY USER
# =========================================================================
parser = argparse.ArgumentParser()
parser.add_argument(
    "--cohort", type=str, choices=["NABEC", "HBCC"], help="Cohort being analyzed", required=True,
)
parser.add_argument(
    "--work_dir", type=str, help="Path to working directory", required=True,
)
parser.add_argument(
    "--chunk_size", type=int, help="Number of promoters per chunk in OLS", required=True,
)
parser.add_argument(
    "--corr_threshold", type=float, help="Filtering cutoff for Pearson's correlation", required=True,
)
parser.add_argument(
    "--target_column", type=str, help="Name of target feature", required=True,
)
parser.add_argument(
    "--vif", type=str, help="Whether to use VIF for filtering", required=True
)
parser.add_argument(
    "--vif_threshold", type=int, help="VIF threshold if running VIF", required=True,
)
args = parser.parse_args()
COHORT = args.cohort
WORK_DIR = args.work_dir
CHUNK_SIZE = args.chunk_size
CORR_THRESHOLD = args.corr_threshold
TARGET_COLUMN = args.target_column
RUN_VIF = True if args.vif.lower() == "true" else False
VIF_THRESHOLD = args.vif_threshold
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# HELPER FUNCTIONS
# =========================================================================
def make_safe_varname(name):
    return re.sub(r'[^0-9a-zA-Z_]', '_', name)

def filter_correlated(df, threshold=CORR_THRESHOLD):
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    corr_matrix = df[numeric_cols].corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    drop_cols = [column for column in upper.columns if any(upper[column] > threshold)]
    return df.drop(columns=drop_cols, errors='ignore')

def calculate_vif(df_subset):
    X = df_subset.select_dtypes(include=[np.number])
    vif_data = pd.DataFrame()
    vif_data["Variable"] = X.columns
    vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    return vif_data


# =========================================================================
# LOAD DATA AND PROCESS CHUNKS
# =========================================================================
for suffix in ["", "_under80"]:
    # Make directories
    OUTPUT_DIR = f"{PREPROCESSING_DIR}/{COHORT}/promoter_data{suffix}"
    CHUNKS_DIR = os.path.join(OUTPUT_DIR, "chunks")
    SIG_DIR = os.path.join(OUTPUT_DIR, "sig_chunks")
    LOG_FILE = os.path.join(SIG_DIR, "pipeline_log.txt")

    # Load data
    INPUT_FILE = f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_covariates_merged{suffix}.txt"
    for i, input_file in enumerate(INPUT_FILE):
        df = pd.read_csv(
            INPUT_FILE, 
            sep="\t",
        )
        covariates = [f"PC{pc_num}" for pc_num in range(1, 21)] + ["sex_binary", "PMI"]

        # Identify gene columns
        non_gene_cols = [TARGET_COLUMN] + covariates
        excluded_cols = set(["sample"] + non_gene_cols)
        gene_cols = [col for col in df.columns if col not in excluded_cols]

    # Process chunks
    total_sig_rows = 0
    num_chunks = 0
    log_entries = []

    for i, start in enumerate(range(0, len(gene_cols), CHUNK_SIZE)):
        chunk_genes = gene_cols[start:start + CHUNK_SIZE]
        chunk_df = df["sample"].to_frame().join(df[non_gene_cols + chunk_genes])

        filtered_chunk_genes_df = filter_correlated(chunk_df[chunk_genes], threshold=CORR_THRESHOLD)
        filtered_genes = filtered_chunk_genes_df.columns.tolist()

        if RUN_VIF:
            vif_df = calculate_vif(filtered_chunk_genes_df)
            final_genes = vif_df[vif_df["VIF"] < VIF_THRESHOLD]["Variable"].tolist()
        else:
            final_genes = filtered_genes

        chunk_results = []

        for gene in final_genes:
            safe_gene = make_safe_varname(gene)
            chunk_df = chunk_df.rename(columns={gene: safe_gene})

            formula = f"{TARGET_COLUMN} ~ {safe_gene} + {' + '.join(covariates)}"
            model = smf.ols(formula, data=chunk_df).fit()

            gene_coef = model.params[safe_gene]
            gene_pval = model.pvalues[safe_gene]
            gene_stderr = model.bse[safe_gene]
            gene_tval = model.tvalues[safe_gene]
            conf_int = model.conf_int().loc[safe_gene]
            conf_low = conf_int[0]
            conf_high = conf_int[1]

            chunk_results.append({
                "Gene": gene,
                "Coefficient": gene_coef,
                "StdErr": gene_stderr,
                "TValue": gene_tval,
                "pval": gene_pval,
                "CI_L95": conf_low,
                "CI_U95": conf_high,
                "R2": model.rsquared,
                "Adj_R2": model.rsquared_adj,
                "AIC": model.aic,
                "BIC": model.bic,
                "F_statistic": model.fvalue,
                "F_pvalue": model.f_pvalue,
                "n_obs": int(model.nobs)
            })

        chunk_results_df = pd.DataFrame(chunk_results)
        chunk_file = os.path.join(CHUNKS_DIR, f"chunk_{i}_results.tsv")
        chunk_results_df.to_csv(
            chunk_file, 
            sep="\t", 
            index=False,
        )

        sig_chunk_df = chunk_results_df[chunk_results_df["pval"] < 0.05]
        sig_chunk_file = os.path.join(SIG_DIR, f"chunk_{i}_sig_results.tsv")
        sig_chunk_df.to_csv(
            sig_chunk_file, 
            sep="\t", 
            index=False,
        )
        sig_count = sig_chunk_df.shape[0]
        total_sig_rows += sig_count
        num_chunks += 1

        log_entries.append(f"Chunk {i}: total_genes={len(chunk_genes)}, filtered_genes={len(filtered_genes)}, final_genes={len(final_genes)}, significant={sig_count}")

        print(f"Chunk {i} done. Genes: {len(chunk_genes)} => Correlated filter: {len(filtered_genes)} => " +
            (f"VIF filter: {len(final_genes)} => " if RUN_VIF else "") +
            f"Saved: {chunk_file} | Significant: {sig_chunk_file}")

    # Write to log file
    with open(LOG_FILE, "w") as log:
        log.write("Regression Pipeline Log\n")
        log.write(f"Total Chunks: {num_chunks}\n")
        log.write(f"Total Significant Results (p<0.05): {total_sig_rows}\n\n")
        for entry in log_entries:
            log.write(entry + "\n")

    print(f"\nLog written to: {LOG_FILE}")
