#!/bin/env python

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import glob
import os
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
# MERGE ALL SIGNIFICANT CHUNKS
# =========================================================================
for suffix in ["", "_under80"]:
    # Define file paths
    INPUT_FILE = f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_covariates_merged{suffix}.txt"
    SIG_DIR = f"{PREPROCESSING_DIR}/{COHORT}/promoter_data{suffix}/sig_chunks"
    OUTPUT_ALL_SIG = os.path.join(SIG_DIR, "all_significant_results.csv")
    OUTPUT_GENE_LIST = os.path.join(SIG_DIR, "significant_promoter_list.txt")

    # Merge chunks
    sig_files = glob.glob(os.path.join(SIG_DIR, "chunk_*_sig_results.tsv"))
    all_dfs = [pd.read_csv(f, sep="\t") for f in sig_files]

    if all_dfs:
        merged_df = pd.concat(all_dfs, ignore_index=True)
        merged_df.to_csv(
            OUTPUT_ALL_SIG, 
            index=False,
        )

        promoter_list = merged_df["Gene"].drop_duplicates().sort_values().tolist()
        sig_promoters = set(promoter_list)
        with open(OUTPUT_GENE_LIST, "w") as out:
            for promoter in sig_promoters:
                out.write(promoter + "\n")

        print(f"✅ Merged {len(sig_files)} files into: {OUTPUT_ALL_SIG}")
        print(f"✅ Saved {len(sig_promoters)} unique promoters to: {OUTPUT_GENE_LIST}")

        # Check if promoters match original data columns
        df_full = pd.read_csv(
            INPUT_FILE, 
            sep="\t",
        )
        excluded_cols = ["sample", "age", "PMI", "sex_binary"]
        excluded_cols += [f"PC{pc_num}" for pc_num in range(1, 21)]
        original_promoters = set(df_full.columns) - set(excluded_cols)

        matched = sig_promoters & original_promoters
        missing = sig_promoters - original_promoters

        print(f"🔍 Match Check:")
        print(f"✅ Matched promoters: {len(matched)}")
        print(f"❌ Missing promoters: {len(missing)}")
        if missing:
            print("Missing promoter names:")
            for promoter in sorted(missing):
                print("  ", promoter)
        print()
        
    else:
        print("⚠️ No significant chunk files found.")
