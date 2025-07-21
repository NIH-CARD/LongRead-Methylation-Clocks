#!/usr/bin/env python3

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import pandas as pd
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
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# COMBINE HAPLOTYPES AT THE SAME CpG
# =========================================================================
# Only execute for HBCC since NABEC data are already unphased.
if COHORT == "HBCC":
    # Load participant-level haplotype data
    path_hap1 = f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_CpGs_hap1.tsv"
    path_hap2 = f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_CpGs_hap2.tsv"
    df_hap1 = pd.read_csv(path_hap1, sep="\t")
    df_hap2 = pd.read_csv(path_hap2, sep="\t")

    # Reformat participant IDs
    df_hap1.columns = [col.replace("GRCh38_1_", "") for col in df_hap1.columns]
    df_hap2.columns = [col.replace("GRCh38_2_", "") for col in df_hap2.columns]

    # Combine haplotypes together by adding read counts and coverage
    df_combined = pd.concat([df_hap1, df_hap2], ignore_index=True)
    df_combined = df_combined.groupby(["#chrom", "start", "end"], as_index=False).sum()

    # Re-compute averages
    avg_col_names = [col for col in df_combined.columns.values if "modFraction" in col]
    reads_col_names = [col.replace("modFraction", "modReads") for col in avg_col_names]
    cov_col_names = [col.replace("modFraction", "validCov") for col in avg_col_names]
    df_combined[avg_col_names] = df_combined[reads_col_names].values / df_combined[cov_col_names].replace(0, 1).values * 100

    # Save results
    df_combined.to_csv(
        f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_CpGs_merged.tsv", 
        sep="\t", 
        index=False,
    )
