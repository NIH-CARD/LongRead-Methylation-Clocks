#!/bin/env python

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
    "--work_dir", type=str, help="Path to working directory", required=True,
)
args = parser.parse_args()
WORK_DIR = args.work_dir


# =========================================================================
# LOAD PROMOTER - GENE MAPPING
# =========================================================================
df_promtoer_gene_mapping = pd.read_csv(
    f"{WORK_DIR}/results/preprocessing/promoter_gene_mapping.txt", 
    sep="\t",
    header=None,
    names=["Promoter", "Gene"],
)


# =========================================================================
# KEEP ONLY SIGNIFICANT GENES
# =========================================================================
significant_genes_nabec = pd.read_csv(
    f"{WORK_DIR}/results/NABEC_significant_genes.txt",
    header=None,
    names=["Gene"],
)
significant_genes_hbcc = pd.read_csv(
    f"{WORK_DIR}/results/HBCC_significant_genes.txt",
    header=None,
    names=["Gene"],
)
significant_genes = pd.concat(significant_genes_nabec, significant_genes_hbcc)
significant_genes.drop_duplicates(inplace=True)
df_promtoer_gene_mapping = df_promtoer_gene_mapping.merge(significant_genes, on="Gene")


# =========================================================================
# LOAD PROMOTER DATA AND KEEP ONLY PROMOTERS FOR SIGNIFICANT GENES
# =========================================================================
df_promoter_positions = pd.read_csv(
    f"{WORK_DIR}/new_promoters_MAY2025/promoter_positions_2kbp.bed", 
    sep="\t", 
    header=None,
)
df_promoter_positions = df_promoter_positions[[0, 1, 2, 3]]
df_promoter_positions.columns = ["chr", "start", "end", "promoter"]
df_promoter_positions = df_promoter_positions.merge(df_promtoer_gene_mapping[["Promoter"]], on="Promoter")


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
