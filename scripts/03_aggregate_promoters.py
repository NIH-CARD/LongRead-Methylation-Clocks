#!/usr/bin/env python3

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import numpy as np
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
parser.add_argument(
    "--data_dir", type=str, help="Path to data directory", required=True,
)
parser.add_argument(
    "--cov_reads", type=int, help="Number of reads for each individual for each CpG to use for QC", required=True,
)
parser.add_argument(
    "--cov_ratio", type=float, help="Fraction of people at or above cov_reads threshold for a given CpG to determine whether that site is used", required=True,
)
args = parser.parse_args()
COHORT = args.cohort
WORK_DIR = args.work_dir
DATA_DIR = args.data_dir
COV_THRESHOLD_READS = args.cov_reads
COV_THRESHOLD_RATIO = args.cov_ratio
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# LOAD PROMOTER DATA
# =========================================================================
df_promoter_positions = pd.read_csv(
    f"{DATA_DIR}/promoter_positions_2kbp.bed", 
    sep="\t", 
    header=None,
)
df_promoter_positions = df_promoter_positions[[0, 1, 2, 3]]
df_promoter_positions.columns = ["chr", "start", "end", "promoter"]


# =========================================================================
# LOAD PARTICIPANT-LEVEL CpG DATA
# =========================================================================
df_cpgs = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_CpGs_merged.tsv", 
    sep="\t", 
    low_memory=False,
)
df_cpgs.rename(columns={"#chrom":"chr"}, inplace=True)


# =========================================================================
# AGGREGATE CpG COVERAGE AND READS ACROSS EACH SPECIFIED PROMOTER REGION
# =========================================================================
promoter_data = []

for _, promoter_row in tqdm(df_promoter_positions.iterrows(), total=len(df_promoter_positions), desc="Aggregating promoters"):
    promoter_chr = promoter_row["chr"]
    promoter_start = promoter_row["start"]
    promoter_end = promoter_row["end"]
    promoter_name = promoter_row["promoter"]

    ## Weighted average after dropping low-coverage sites
    # Filter to matching only include CpGs within the current promoter
    matching_data = df_cpgs[
        (df_cpgs["chr"] == promoter_chr) &
        (df_cpgs["start"] >= promoter_start) &
        (df_cpgs["start"] <= promoter_end)
    ]
    
    ## Aggregate reads and coverage data
    if not matching_data.empty:
        # Keep reads and coverage columns, remove the rest
        cols_to_drop = list(df_cpgs.columns[4::3].values) + ["chr", "start", "end"]
        promoter = matching_data.drop(columns=cols_to_drop)
        promoter = promoter.apply(pd.to_numeric, errors="coerce").fillna(0)

        # Remove rows determined to have low coverage based on % of participants whose coverage is below a defined threshold
        cov_cols = promoter.columns[::2]
        below_filter = promoter[cov_cols] < COV_THRESHOLD_READS
        percent_below_filter = below_filter.sum(axis=1) / len(cov_cols)
        promoter = promoter[percent_below_filter < 0.05]

        if len(promoter) > 0:
            # Sum across columns
            promoter = promoter.sum()
            promoter = promoter.to_frame().T

            # Create a new DataFrame with the averages, ensuring there are no entries with denominator of 0
            promoter.iloc[:, ::2] = promoter.iloc[:, ::2].replace(0, 1)
            arr_avg = promoter.iloc[:, 1::2].values / promoter.iloc[:, 0::2].values * 100

            # Reformat promoter data
            colnames = ["_".join(col.split("_")[:-1]).replace('-', '_') for col in promoter.columns[::2]]
            promoter = pd.DataFrame(arr_avg, columns=colnames, index=[promoter_name])
            promoter_data.append(promoter)

df_promoter = pd.concat(promoter_data)
df_promoter = df_promoter.T.rename_axis('sample').rename_axis(columns=None).reset_index()


# =========================================================================
# SAVE MERGED PROMOTER DATA
# =========================================================================
with open(f'{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoter_columns.txt', 'w') as file:
    for column in df_promoter.columns:
        file.write(f"{column}\n")

df_promoter.to_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_merged_2kb.tsv",
    sep="\t", 
    index=False,
)
