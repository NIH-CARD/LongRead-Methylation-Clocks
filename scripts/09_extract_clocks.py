#!/usr/bin/env python3

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import numpy as np
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
parser.add_argument(
    "--data_dir", type=str, help="Path to data directory", required=True,
)
parser.add_argument(
    "--promoter_length", type=int, help="Length of base pair window for each promoter", required=True,
)
args = parser.parse_args()
COHORT = args.cohort
WORK_DIR = args.work_dir
DATA_DIR = args.data_dir
PROMOTER_LENGTH = args.promoter_length
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# DEFINE HELPER FUNCTIONS
# =========================================================================
def find_closest_promoter(df_clock, df_promoter_ranges, df_promoters):
    """ Find the closest promoter to each clock site """
    log_entries = []
    subset_promoters = []

    for _, clock in df_clock.iterrows():
        # Filter promoters that are on the same chromosome
        same_chrom_promoters = df_promoter_ranges[df_promoter_ranges["chrom"] == clock["chrom"]]
        same_chrom_promoters = same_chrom_promoters[same_chrom_promoters["promoter"].isin(df_promoters.columns.values)]

        # Find the nearest promoter by computing the distance
        same_chrom_promoters = same_chrom_promoters.assign(
            distance=same_chrom_promoters.apply(
                lambda row: abs(row["center"] - clock["start"]), axis=1
            )
        )
        nearest_promoter = same_chrom_promoters.loc[same_chrom_promoters["distance"].idxmin()]
        distance = nearest_promoter["distance"] - PROMOTER_LENGTH / 2
        if distance <= 0:
            distance = 0
            relation = "within"
        else:
            relation = "next_to"

        # Save log entry
        log_entries.append({
            "CGid": clock["CGid"],
            "chrom": clock["chrom"],
            "clock_start": clock["start"],
            "promoter_name": nearest_promoter["promoter"],
            "promoter_start": nearest_promoter["start"],
            "promoter_end": nearest_promoter["end"],
            "relation": relation,
            "distance": distance,
        })

        # Save subset of promoters
        subset_promoters.append(nearest_promoter["promoter"])

    # Create log dataframe and get set of unique promoters
    df_log = pd.DataFrame(log_entries)
    subset_promoters = set(subset_promoters)

    return df_log, subset_promoters

def find_closest_cpgs(row, df_full):
    same_chr = df_full[df_full["chrom"] == row["chrom"]]
    if same_chr.empty:
        return None
    distance = (same_chr["start"] - row["start"]).abs()
    return same_chr.loc[distance.idxmin()]


# =========================================================================
# FIND PROMOTERS OVERLAPPING WITH CLOCK SITES
# =========================================================================
# Pull individual level data for unique promoters
df_promoters = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_merged_2kb.tsv", 
    sep="\t",
)
included_promoters = df_promoters.columns.values[1:]

# Pull promoter positions, only keeping those that are in the current cohort's promoter list
df_promoter_ranges = pd.read_csv(
    f'{DATA_DIR}/promoter_positions_2kbp.bed', 
    sep="\t", 
    header=None,
)
df_promoter_ranges.columns = ["chrom", "start", "end", "promoter", "?", "?"]
df_promoter_ranges = df_promoter_ranges[["promoter", "chrom", "start", "end"]]
df_promoter_ranges["center"] = (df_promoter_ranges["start"] + df_promoter_ranges["end"]) / 2
df_promoter_ranges = df_promoter_ranges[df_promoter_ranges["promoter"].isin(included_promoters)].reset_index(drop=True)

for clock_name in ["clock1", "clock2", "clock3"]:
    df_clock = pd.read_csv(
        f"{WORK_DIR}/clock_positions/{clock_name}_CGid_positions.bed", 
        sep="\t", 
        header=None,
    )
    df_clock.columns = ["chrom", "start", "end", "?", "CGid", "?"]
    df_clock = df_clock[["chrom", "start", "end", "CGid"]]
    df_log, subset_promoters = find_closest_promoter(df_clock, df_promoter_ranges, df_promoters)

    df_promoters_subset = df_promoters[["sample"] + list(subset_promoters)]

    # Save data
    df_log.to_csv(
        f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_{clock_name.upper()}_PROMOTER_LOG.txt", 
        index=False, 
        sep="\t",
    )
    df_promoters_subset.to_csv(
        f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_{clock_name.upper()}_PROMOTER_SUBSET.txt", 
        index=False, 
        sep="\t",
    )


# =========================================================================
# FIND CpGs OVERLAPPING WITH CLOCK SITES
# =========================================================================
for clock_name in ["clock1", "clock2", "clock3"]:
    # Read clock positions
    df_clock_positions = pd.read_csv(
        f"{WORK_DIR}/clock_positions/{clock_name}_CGid_positions.bed", 
        sep="\t", 
        header=None,
    )
    df_clock_positions.rename(columns={0:"chrom", 1:"start", 2:"end"}, inplace=True)
    df_clock_positions = df_clock_positions[["chrom", "start", "end"]]

    # Read haplotype data
    df_hap1 = pd.read_csv(
        f"{WORK_DIR}/clock_overlap/{COHORT}_{clock_name}_hap1.tsv", 
        sep="\t",
        low_memory=False,
    ).rename(columns={"#chrom":"chrom"})
    df_hap2 = pd.read_csv(
        f"{WORK_DIR}/clock_overlap/{COHORT}_{clock_name}_hap2.tsv", 
        sep="\t",
        low_memory=False,
    ).rename(columns={"#chrom":"chrom"})

    # Remove duplicate rows
    df_hap1 = df_hap1.drop_duplicates(keep="first")
    df_hap2 = df_hap2.drop_duplicates(keep="first")

    # Reformat participant IDs
    df_hap1.columns = [col.replace("GRCh38_1_", "") for col in df_hap1.columns]
    df_hap2.columns = [col.replace("GRCh38_2_", "") for col in df_hap2.columns]

    # Combine haplotypes together by adding read counts and coverage
    df_cpgs = pd.concat([df_hap1, df_hap2], ignore_index=True)
    df_cpgs = df_cpgs.groupby(["chrom", "start", "end"], as_index=False).sum()

    # Re-compute averages
    avg_col_names = [col for col in df_cpgs.columns.values if "modFraction" in col]
    reads_col_names = [col.replace("modFraction", "modReads") for col in avg_col_names]
    cov_col_names = [col.replace("modFraction", "validCov") for col in avg_col_names]
    df_cpgs[avg_col_names] = df_cpgs[reads_col_names].values / df_cpgs[cov_col_names].replace(0, 1).values * 100

    # Filter to only include clock CpGs
    df_cpgs = df_clock_positions.apply(find_closest_cpgs, axis=1, df_full=df_cpgs)
    df_cpgs = df_cpgs.copy()
    df_cpgs["index"] = df_cpgs["chrom"] + ":" + df_cpgs["start"].astype(str)
    df_cpgs = df_cpgs.set_index("index")
    df_cpgs.index.name = None
    df_cpgs.drop(columns=["chrom", "start", "end"], inplace=True)
    df_cpgs = df_cpgs.loc[:, df_cpgs.columns.str.contains('modFraction')]
    df_cpgs.columns = df_cpgs.columns.str.replace("_modFraction", "").str.replace("-", "_")
    df_cpgs = df_cpgs.T.reset_index().rename(columns={"index": "ID"})
    df_cpgs.to_csv(
        f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_{clock_name.upper()}_CpG_SUBSET.tsv",
        index=False,
        sep="\t",
    )
