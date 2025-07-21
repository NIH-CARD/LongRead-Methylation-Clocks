#!/bin/env python

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import shap
import statsmodels.formula.api as smf
from tqdm import tqdm


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
args = parser.parse_args()
COHORT = args.cohort
WORK_DIR = args.work_dir
RESULTS_DIR = f"{WORK_DIR}/results"
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# LOAD MODEL AND PARTICIPANT-LEVEL DATA
# =========================================================================
model = joblib.load(f"{RESULTS_DIR}/{COHORT}_promoter/model.joblib")
X_train = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_train_promoter.tsv", 
    sep="\t",
)
y = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_pheno.tsv", 
    sep="\t",
)
df_train = y.merge(X_train, on="ID")
df_train.drop(columns=["ID", "PMI", "sex_binary"], inplace=True)
X_train.drop(columns="ID", inplace=True)


# =========================================================================
# LOAD PROMOTER - GENE MAPPING
# =========================================================================
df_promtoer_gene_mapping = pd.read_csv(
    f"{PREPROCESSING_DIR}/promoter_gene_mapping.txt", 
    sep="\t",
    header=None,
    names=["Promoter", "Gene"],
)


# =========================================================================
# FIND PROMOTERS THAT WERE INCLUDED IN THE MODEL WITH NONZERO COEFFICIENTS
# =========================================================================
promoter_names = X_train.columns.values[2:]
promoter_coefs = model.coef_[2:]
df_coefs = pd.DataFrame({
    "Promoter": promoter_names,
    "Coefficient": promoter_coefs,
})
df_coefs = df_coefs[df_coefs.Coefficient != 0]
relevant_promoters = df_coefs.Promoter

df_gene = df_coefs.merge(df_promtoer_gene_mapping, on="Promoter")
unique_genes = list(set(df_gene["Gene"].tolist()))
with open(f"{RESULTS_DIR}/{COHORT}_full_genes.txt", "w") as file:
    for gene in unique_genes:
        file.write(f"{gene}\n")


# =========================================================================
# RUN REGRESSIONS TO FIND SIGNIFICANCE OF RELATIONSHIP WITH AGE
# =========================================================================
df_results = []
for promoter_name in tqdm(relevant_promoters, total=len(relevant_promoters), desc="Running regressions for each promoter"):
    clean_promoter_name = re.sub(r'[^0-9a-zA-Z_]', '_', promoter_name)
    df_train.rename(columns={promoter_name: clean_promoter_name}, inplace=True)

    formula = f"PHENO ~ {clean_promoter_name}"
    promoter_model = smf.ols(formula, data=df_train).fit()

    coef = promoter_model.params[clean_promoter_name]
    stderr = promoter_model.bse[clean_promoter_name]
    pval = promoter_model.pvalues[clean_promoter_name]
    conf_int = promoter_model.conf_int().loc[clean_promoter_name]
    conf_low = conf_int[0]
    conf_high = conf_int[1]
    df_results.append(pd.DataFrame({
        "Promoter": [promoter_name],
        "Coefficient": [coef],
        "StdErr": [stderr],
        "pval": [pval],
        "CI_L95": [conf_low],
        "CI_U95": [conf_high],
        "R2": [promoter_model.rsquared],
        "Adj_R2": [promoter_model.rsquared_adj],
        "AIC": [promoter_model.aic],
        "BIC": [promoter_model.bic],
        "F_statistic": [promoter_model.fvalue],
        "F_pvalue": [promoter_model.f_pvalue],
        "n_obs": [int(promoter_model.nobs)],
    }))
df_results = pd.concat(df_results)
df_results.sort_values(by="pval", inplace=True)
df_results["Bonferroni_threshold"] = 0.05 / len(df_results)
df_results = df_results[df_results["pval"] < df_results["Bonferroni_threshold"]]

df_significant_gene = df_results.merge(df_promtoer_gene_mapping, on="Promoter")
unique_genes = list(set(df_significant_gene["Gene"].tolist()))
with open(f"{RESULTS_DIR}/{COHORT}_significant_genes.txt", "w") as file:
    for gene in unique_genes:
        file.write(f"{gene}\n")
