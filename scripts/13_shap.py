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
    "--cohort", type=str, choices=["NABEC", "HBCC"], help="Cohort being analyzed", required=True,
)
parser.add_argument(
    "--work_dir", type=str, help="Path to working directory", required=True,
)
parser.add_argument(
    "--n_shap", type=int, help="Number of SHAP values to include", required=True,
)
args = parser.parse_args()
COHORT = args.cohort
WORK_DIR = args.work_dir
N_SHAP = args.n_shap
RESULTS_DIR = f"{WORK_DIR}/results"
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"


# =========================================================================
# LOAD MODEL AND PARTICIPANT-LEVEL DATA
# =========================================================================
model = joblib.load(f"{RESULTS_DIR}/{COHORT}_promoter/model.joblib")
X_test = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_test_promoter.tsv", 
    sep="\t",
)
y = pd.read_csv(
    f"{PREPROCESSING_DIR}/{COHORT}/{COHORT}_pheno.tsv", 
    sep="\t",
)
df_test = y.merge(X_test, on="ID")
df_test.drop(columns=["ID", "PMI", "sex_binary"], inplace=True)
X_test.drop(columns="ID", inplace=True)


# =========================================================================
# CALCULATE SHAP VALUES
# =========================================================================
# All SHAP values
explainer = shap.Explainer(model, X_test)
shap_values = explainer(X_test)

# Subset top N_SHAP values for plots where the "Sum of other features" would dilute visualization
mean_abs_shap = np.abs(shap_values.values).mean(axis=0)
top_idx = np.argsort(mean_abs_shap)[::-1][:N_SHAP]
shap_top = shap.Explanation(
    values=shap_values.values[:, top_idx],
    base_values=shap_values.base_values,
    data=shap_values.data[:, top_idx],
    feature_names=[shap_values.feature_names[i] for i in top_idx]
)


# =========================================================================
# SAVE SHAP PLOTS
# =========================================================================
shap.plots.beeswarm(shap_top, max_display=N_SHAP, show=False)
plt.savefig(f"{RESULTS_DIR}/{COHORT}_promoter/shap_{COHORT}_beeswarm.png", dpi=600, bbox_inches="tight")
plt.figure()

shap.plots.bar(shap_values, max_display=N_SHAP+1, show=False)
plt.savefig(f"{RESULTS_DIR}/{COHORT}_promoter/shap_{COHORT}_bar.png", dpi=600, bbox_inches="tight")
plt.figure()


# # =========================================================================
# # COMBINE BEESWARM AND BAR PLOTS
# # =========================================================================
# whitespace_pixels = 250
# colorbar_width = 600
# beeswarm_feature_label_width = 1000

# # Load beeswarm and bar plots
# arr_beeswarm = np.array(Image.open(f"{RESULTS_DIR}/{COHORT}_promoter/shap_{COHORT}_beeswarm.png"))
# arr_barplot = np.array(Image.open(f"{RESULTS_DIR}/{COHORT}_promoter/shap_{COHORT}_bar.png"))

# # Resize beeswarm plot, keeping aspect ratio constant
# beeswarm_dim0 = arr_barplot.shape[0]
# beeswarm_dim1 = int(arr_beeswarm.shape[1] * (arr_barplot.shape[0] / arr_beeswarm.shape[0]))
# img_beeswarm = Image.fromarray(arr_beeswarm)
# img_beeswarm = img_beeswarm.resize((beeswarm_dim1, beeswarm_dim0), Image.LANCZOS)
# arr_beeswarm = np.array(img_beeswarm)

# # Add whitespace at the bottom of the beeswarm plot to account for extra "sum of features" row in bar plot
# white_rows = np.ones((whitespace_pixels, arr_beeswarm.shape[1], arr_beeswarm.shape[2]), dtype=np.uint8) * 255
# arr_beeswarm = np.vstack((arr_beeswarm, white_rows))

# # Resize bar plot, keeping aspect ratio the same, to match whitespace added in beeswarm plot
# barplot_dim0 = arr_barplot.shape[0] + whitespace_pixels
# barplot_dim1 = arr_barplot.shape[1] + int(whitespace_pixels*(arr_barplot.shape[1] / arr_barplot.shape[0]))
# img_barplot = Image.fromarray(arr_barplot)
# img_barplot = img_barplot.resize((barplot_dim1, barplot_dim0), Image.LANCZOS)
# arr_barplot = np.array(img_barplot)

# # Stitch together beeswarm colorbar, beeswarm plot, and bar plot
# beeswarm_colorbar = arr_beeswarm[:, -colorbar_width:]
# beeswarm_plot = arr_beeswarm[:, beeswarm_feature_label_width:-colorbar_width]
# arr_merged = np.concatenate((beeswarm_colorbar, beeswarm_plot, arr_barplot), axis=1)

# # Save merged image
# Image.fromarray(arr_merged).save(
#     f"{RESULTS_DIR}/{COHORT}_promoter/shap_{COHORT}_merged.png",
#     dpi=(600, 600),
# )
