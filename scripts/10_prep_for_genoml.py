#!/bin/env python

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import pandas as pd
from sklearn.model_selection import train_test_split


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
# MERGE PROMOTER DATA AND METADATA AND PERFORM TRAIN/TEST SPLITS
# =========================================================================
for suffix in ["", "_under80"]:
    # Load and merge promoter data and metadata
    sig_promoters_path = f"{PREPROCESSING_DIR}/{COHORT}/promoter_data{suffix}/sig_chunks/significant_promoter_list.txt"
    with open(sig_promoters_path, 'r') as file:
        sig_promoters_list = [line.strip() for line in file.readlines()]
    df_promoter = pd.read_csv(
        f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_merged_2kb.tsv", 
        sep="\t",
    )
    df_promoter = df_promoter[["sample"] + sig_promoters_list]
    df_promoter.rename(columns={"sample":"ID"}, inplace=True)

    df_metadata = pd.read_csv(f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_metadata.tsv", sep="\t")
    df_metadata.rename(columns={"sample":"ID", "age":"PHENO"}, inplace=True)
    df_metadata = df_metadata[["ID", "PMI", "sex_binary", "PHENO"]]

    df_merged = df_metadata.merge(df_promoter, on="ID")
    if suffix == "_under80":
        df_merged = df_merged[df_merged["PHENO"] <= 80]

    # Save pheno file
    if suffix == "":
        df_pheno = df_metadata[["ID", "PHENO"]]
        df_pheno.to_csv(
            f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_pheno.tsv", 
            sep="\t", 
            index=False,
        )

    # Train/test split
    y = df_merged[["PHENO"]]
    X = df_merged.drop(columns=["PHENO"])
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=3)
    X_train.to_csv(
        f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_train_promoter{suffix}.tsv", 
        sep="\t", 
        index=False,
    )
    X_test.to_csv(
        f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_test_promoter{suffix}.tsv", 
        sep="\t", 
        index=False,
    )

    # Apply train/test split on clock ocerlap datasets
    if suffix == "":
        ids_train = X_train["ID"].tolist()
        ids_test = X_test["ID"].tolist()
        for clock_name in ["clock1", "clock2", "clock3"]:
            # Promoter overlap
            df_clock_promoter = pd.read_csv(
                f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_{clock_name.upper()}_PROMOTER_SUBSET.txt",
                sep="\t",
            )
            df_clock_promoter.rename(columns={"sample":"ID"}, inplace=True)
            df_clock_promoter_train = df_clock_promoter[df_clock_promoter["ID"].isin(ids_train)]
            df_clock_promoter_test = df_clock_promoter[df_clock_promoter["ID"].isin(ids_test)]
            df_clock_promoter_train.to_csv(
                f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_train_promoter_{clock_name}.tsv", 
                sep="\t", 
                index=False,
            )
            df_clock_promoter_test.to_csv(
                f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_test_promoter_{clock_name}.tsv", 
                sep="\t", 
                index=False,
            )

            # CpG overlap
            df_clock_cpgs = pd.read_csv(
                f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_{clock_name.upper()}_CpG_SUBSET.tsv",
                sep="\t",
            )
            df_clock_cpgs_train = df_clock_cpgs[df_clock_cpgs["ID"].isin(ids_train)]
            df_clock_cpgs_test = df_clock_cpgs[df_clock_cpgs["ID"].isin(ids_test)]
            df_clock_cpgs_train.to_csv(
                f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_train_CpG_{clock_name}.tsv", 
                sep="\t", 
                index=False,
            )
            df_clock_cpgs_test.to_csv(
                f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_test_CpG_{clock_name}.tsv", 
                sep="\t", 
                index=False,
            )
