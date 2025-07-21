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
    "--work_dir", type=str, help="Path to working directory", required=True,
)
args = parser.parse_args()
WORK_DIR = args.work_dir
RESULTS_DIR = f"{WORK_DIR}/results"
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


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
dfs = {}
for cohort in ["NABEC", "HBCC"]:
    for suffix in ["_promoter", "_promoter_under80"]:
        model = joblib.load(f"{RESULTS_DIR}/{cohort + suffix}/model.joblib")
        promoter_names = pd.read_csv(
            f"{PREPROCESSING_DIR}/{cohort}/{cohort}_train{suffix}.tsv", 
            sep="\t",
        ).columns.values[3:]
        promoter_coefs = model.coef_[2:]
        df_coefs = pd.DataFrame({
            "Promoter": promoter_names,
            "Coefficient": promoter_coefs,
        })
        df_coefs = df_coefs[df_coefs.Coefficient != 0]
        relevant_promoters = df_coefs.Promoter
        df_gene = df_coefs.merge(df_promtoer_gene_mapping, on="Promoter")
        dfs[cohort + suffix] = df_gene
df_nabec = dfs["NABEC_promoter"]
df_nabec_under80 = dfs["NABEC_promoter_under80"]
df_hbcc = dfs["HBCC_promoter"]
df_hbcc_under80 = dfs["HBCC_promoter_under80"]


# =========================================================================
# PRINT OUT OVERLAPS
# =========================================================================
print(f"Promoters in NABEC:                             {len(df_nabec)}")
print(f"Promoters in HBCC:                              {len(df_hbcc)}")
print(f"Promoters in NABEC AND HBCC:                    {len(df_nabec[df_nabec["Promoter"].isin(df_hbcc["Promoter"])])}")
print(f"Promoters in NABEC UNDER 80:                    {len(df_nabec_under80)}")
print(f"Promoters in HBCC UNDER 80:                     {len(df_hbcc_under80)}")
print(f"Promoters in NABEC UNDER 80 AND HBCC UNDER 80:  {len(df_nabec_under80[df_nabec_under80["Promoter"].isin(df_hbcc_under80["Promoter"])])}")
print(f"Promoters in NABEC AND NABEC UNDER 80:          {len(df_nabec[df_nabec["Promoter"].isin(df_nabec_under80["Promoter"])])}")
print(f"Promoters in HBCC AND HBCC UNDER 80:            {len(df_hbcc[df_hbcc["Promoter"].isin(df_hbcc_under80["Promoter"])])}")
print()
print()
print(f"Unique genes in NABEC:                              {len(set(df_nabec["Gene"]))}")
print(f"Unique genes in HBCC:                               {len(set(df_hbcc["Gene"]))}")
print(f"Unique genes in NABEC AND HBCC:                     {len(set(df_nabec.merge(df_hbcc[["Promoter"]], on="Promoter")["Gene"]))}")
print(f"Unique genes in NABEC UNDER 80:                     {len(set(df_nabec_under80["Gene"]))}")
print(f"Unique genes in HBCC UNDER 80:                      {len(set(df_hbcc_under80["Gene"]))}")
print(f"Unique genes in NABEC UNDER 80 AND HBCC UNDER 80:   {len(set(df_nabec_under80.merge(df_hbcc_under80[["Promoter"]], on="Promoter")["Gene"]))}")
print(f"Unique genes in NABEC AND NABEC UNDER 80:           {len(set(df_nabec.merge(df_nabec_under80[["Promoter"]], on="Promoter")["Gene"]))}")
print(f"Unique genes in HBCC AND HBCC UNDER 80:             {len(set(df_hbcc.merge(df_hbcc_under80[["Promoter"]], on="Promoter")["Gene"]))}")


