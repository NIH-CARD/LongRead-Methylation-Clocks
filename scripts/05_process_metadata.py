#!/usr/bin/env python3

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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
# READ AND REFORMAT METADATA
# =========================================================================
df_meta = pd.read_csv(
    f"{WORK_DIR}/meta_data/{COHORT}_meta.tsv", 
    sep="\t",
)
df_meta.columns = df_meta.columns.str.replace("_", "")
if COHORT=="NABEC":
    df_meta = df_meta[df_meta["ancestry"] == "EUR"]
elif COHORT=="HBCC":
    df_meta = df_meta[df_meta["ancestry"].isin(["AFR", "AAC"])]
df_meta["sex_binary"] = np.where(df_meta["sex"].str.lower() == "male", 1, 0)
if len(set(df_meta['group'].tolist())) > 1:
    group_dummies = pd.get_dummies(df_meta['group'], prefix='group').astype(int)
    df_meta = pd.concat([df_meta, group_dummies], axis=1)
df_meta = df_meta.drop(['group'], axis=1)


# =========================================================================
# BIN AGE IN 10-YEAR INCREMENTS
# =========================================================================
bins = list(range(0, 101, 10))
df_meta["AGE_BIN"] = pd.cut(
    df_meta["age"],
    bins=bins,
    labels=[f"{i}–{i+9}" for i in range(0, 100, 10)],
    include_lowest=True,
    right=False           # intervals are left-inclusive, right-exclusive [0, 10)
)
age_counts = df_meta.groupby(["AGE_BIN", "sex"])["sample"].count()
age_counts_pivot = age_counts.unstack("sex").fillna(0)


# =========================================================================
# PLOT AGE FOR EACH SEX
# =========================================================================
fig, ax = plt.subplots(figsize=(10, 6))
age_counts_pivot.plot(kind="bar", ax=ax, color=["#B359BF", "#4D26BF"])
for container in ax.containers:
    ax.bar_label(container, label_type="edge", fontsize=8)
ax.set_xlabel("Age range (years)")
ax.set_ylabel("Count")
ax.set_title(f"{COHORT} Age Distribution by Sex")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Sex")
plt.tight_layout()
plt.savefig(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_Age_Distribution.png", 
    dpi=600,
)


# =========================================================================
# GENERATE CORRELATION MATRIX
# =========================================================================
df_corr = df_meta.drop(columns = [df_meta.columns[0]] + ["sex", "ancestry", "AGE_BIN"])
corr_matrix = df_corr.corr()
mask = np.zeros_like(corr_matrix, dtype=bool)
mask[np.triu_indices_from(mask)] = True


# =========================================================================
# CREATE HEATMAP WITH CORRECT ASPECT RATIO USING A CUSTOM DIVERGING COLORMAP
# =========================================================================
f, ax = plt.subplots(figsize=(11, 9))
cmap = sns.diverging_palette(220, 10, as_cmap=True)
sns.heatmap(
    corr_matrix, 
    mask=mask, 
    cmap=cmap, 
    vmax=1, 
    vmin=-1, 
    center=0,
    square=True, 
    linewidths=1, 
    cbar_kws={"shrink": .5},
    annot=True, 
    fmt=".1f", 
    annot_kws={"size": 6},
)
plt.savefig(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_Metadata_Correlation.png", 
    dpi=600,
)
plt.title(f"Correlation Matrix of Metadata in {COHORT}")


# =========================================================================
# SAVE METADATA
# =========================================================================
columns_to_drop = [col for col in df_meta.columns if col.startswith('group')] + ["sex", "ancestry", "AGE_BIN"]
df_meta.drop(columns=columns_to_drop, inplace=True)
df_meta.to_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_metadata.tsv", 
    sep="\t", 
    index=False,
)
