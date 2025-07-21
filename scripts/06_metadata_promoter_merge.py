#!/usr/bin/env python3

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import pandas as pd


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
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# REPLACE "-" WITH "_" IN THE SAMPLE COLUMN FOR BOTH DATAFRAMES
# =========================================================================
df_promoter = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_merged_2kb.tsv", 
    sep="\t",
)
df_meta = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_metadata.tsv", 
    sep="\t",
)
df_merged = df_meta.merge(df_promoter, on="sample")


# =========================================================================
# SAVE MERGED DATA 
# =========================================================================
df_merged.to_csv(
    f'{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_covariates_merged.txt', 
    index=False, 
    sep="\t",
)
df_merged[df_merged["age"] <= 80].to_csv(
    f'{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_covariates_merged_under80.txt', 
    index=False, 
    sep="\t",
)
with open(f'{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoter_meta_merged_columns.txt', 'w') as file:
    for column in df_merged.columns.tolist():
        file.write(f"{column}\n")


# =========================================================================
# PRINT INFO ON NUMBER OF SAMPLES IN PROMOTER DATA AND METADATA
# =========================================================================
promoter_no_meta = df_promoter[~df_promoter['sample'].isin(df_meta['sample'])]['sample'].tolist()
meta_no_promoter = df_meta[~df_meta['sample'].isin(df_promoter['sample'])]['sample'].tolist()
print(COHORT)
print(f"Number of merged samples:                      {df_merged['sample'].nunique()}")
print(f"Total samples in filtered metadata:            {df_meta['sample'].nunique()}")
print(f"Total samples in promoter data:                {df_promoter['sample'].nunique()}")
print(f"\nSamples in promoter data but not metadata")
print(promoter_no_meta)
print(f"\nSamples in metadata but not promoter data")
print(meta_no_promoter)
print(f"Number of merged samples below 80:             {df_merged[df_merged['age'] <= 80]['sample'].nunique()}")
print()
print()
