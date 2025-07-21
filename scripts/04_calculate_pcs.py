#!/usr/bin/env python3

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


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
    "--n_pcs", type=int, help="Number of PCs to keep", required=True,
)
args = parser.parse_args()
COHORT = args.cohort
WORK_DIR = args.work_dir
N_PCS = args.n_pcs
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# READ PROMOTER DATA -- ROWS AS PARTICIPANTS, COLUMNS AS PROMOTERS
# =========================================================================
df = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_promoters_merged_2kb.tsv", 
    sep="\t",
)
df = df.set_index('sample')


# =========================================================================
# FIT PCA, KEEPING NUMBER OF PCs DEFINED BY USER
# =========================================================================
pca = PCA(n_components=N_PCS, random_state=3)
arr_pcs = pca.fit_transform(df)


# =========================================================================
# SAVE PCs FOR EACH PARTICIPANT
# =========================================================================
df_pcs = pd.DataFrame(
    arr_pcs, 
    index=df.index, 
    columns=[f"methPC{i}" for i in range(1, N_PCS + 1)],
)
df_pcs.reset_index(names="ID", inplace=True)
df_pcs.to_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_pcs.txt", 
    sep="\t", 
    index=False,
)


# =========================================================================
# GENERATE SCREE PLOT
# =========================================================================
variance_explained = pca.explained_variance_ratio_ * 100
plot_data = pd.DataFrame({
    'PC': np.arange(1, N_PCS + 1),
    'VarianceExplained': variance_explained
})
sns.set(rc={'figure.figsize':(12, 8)})
scree_plot = sns.pointplot(
    x='PC', 
    y='VarianceExplained', 
    data=plot_data,
)
sns.set_style("whitegrid")
sns.set_context("talk")
plt.grid(False)
scree_plot.set_title(f'Scree Plot: {COHORT.upper()} PCA')
scree_plot.set_ylabel('Percent (%) Variance Explained')
scree_plot.set_xlabel('Principal Components')
plt.ylim(0, variance_explained.max() + 5)
plt.savefig(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_screeplot.png", 
    dpi=600,
)
plt.clf()
