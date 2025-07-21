#!/bin/env python

# =========================================================================
# IMPORT PACKAGES
# =========================================================================
import argparse
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from PIL import Image, ImageDraw, ImageFont


# =========================================================================
# PARSE PARAMS PASSED BY USER
# =========================================================================
parser = argparse.ArgumentParser()
parser.add_argument(
    "--work_dir", type=str, help="Path to working directory", required=True,
)
args = parser.parse_args()
WORK_DIR = args.work_dir
RESULTS_DIR = f"{WORK_DIR}/results"
PREPROCESSING_DIR = f"{WORK_DIR}/results/preprocessing"
COHORTS = ["HBCC", "NABEC"]


# =========================================================================
# LOAD FONT FOR PLOTS
# =========================================================================
subprocess.run(f"wget -O {RESULTS_DIR}/roboto.zip https://dl.dafont.com/dl/?f=roboto", check=True, shell=True)
subprocess.run(f"unzip -o {RESULTS_DIR}/roboto.zip -d {RESULTS_DIR}/roboto", check=True, shell=True)


# =========================================================================
# FIGURE 4 -- COMBINE BEESWARM AND BAR PLOTS
# =========================================================================
whitespace_pixels = 250
colorbar_width = 600
beeswarm_feature_label_width = 1000
label_padding = 100
font_size = 300

imgs_labeled = []
total_height = 0
first_image_width = None
for i, cohort in enumerate(COHORTS):
    # Load beeswarm and bar plots
    arr_beeswarm = np.array(Image.open(f"{RESULTS_DIR}/{cohort}_promoter/shap_{cohort}_beeswarm.png"))
    arr_barplot = np.array(Image.open(f"{RESULTS_DIR}/{cohort}_promoter/shap_{cohort}_bar.png"))

    # Resize beeswarm plot, keeping aspect ratio constant
    beeswarm_dim0 = arr_barplot.shape[0]
    beeswarm_dim1 = int(arr_beeswarm.shape[1] * (arr_barplot.shape[0] / arr_beeswarm.shape[0]))
    img_beeswarm = Image.fromarray(arr_beeswarm)
    img_beeswarm = img_beeswarm.resize((beeswarm_dim1, beeswarm_dim0), Image.LANCZOS)
    arr_beeswarm = np.array(img_beeswarm)

    # Add whitespace at the bottom of the beeswarm plot to account for extra "sum of features" row in bar plot
    white_rows = np.ones((whitespace_pixels, arr_beeswarm.shape[1], arr_beeswarm.shape[2]), dtype=np.uint8) * 255
    arr_beeswarm = np.vstack((arr_beeswarm, white_rows))

    # Resize bar plot, keeping aspect ratio the same, to match whitespace added in beeswarm plot
    barplot_dim0 = arr_barplot.shape[0] + whitespace_pixels
    barplot_dim1 = arr_barplot.shape[1] + int(whitespace_pixels*(arr_barplot.shape[1] / arr_barplot.shape[0]))
    img_barplot = Image.fromarray(arr_barplot)
    img_barplot = img_barplot.resize((barplot_dim1, barplot_dim0), Image.LANCZOS)
    arr_barplot = np.array(img_barplot)

    # Stitch together beeswarm colorbar, beeswarm plot, and bar plot
    beeswarm_colorbar = arr_beeswarm[:, -colorbar_width:]
    beeswarm_plot = arr_beeswarm[:, beeswarm_feature_label_width:-colorbar_width]
    arr_merged = np.concatenate((beeswarm_colorbar, beeswarm_plot, arr_barplot), axis=1)
    img_merged = Image.fromarray(arr_merged)

    # Find height in pixels for cohort label
    font = ImageFont.truetype(f"{RESULTS_DIR}/roboto/Roboto-Regular.ttf", size=font_size)
    dummy_img = Image.new("RGB", (1, 1))
    draw = ImageDraw.Draw(dummy_img)
    bbox = draw.textbbox((0, 0), cohort, font=font)

    label_height = bbox[3] - bbox[1] + label_padding * 2
    label_width = bbox[2] - bbox[0]

    # Create new image with both the figure and the label
    img_labeled = Image.new("RGBA", (img_merged.width, label_height + img_merged.height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img_labeled)
    draw.text((img_merged.width // 2 - label_width // 2, label_padding), cohort, font=font, fill=(0, 0, 0, 255))
    img_labeled.paste(img_merged, (0, label_height))

    if i == 0:
        first_image_width = img_merged.width
    else:
        new_height = (label_height + img_merged.height) * (first_image_width / img_merged.width)
        img_labeled = img_labeled.resize((first_image_width, int(new_height)), Image.LANCZOS)

    total_height += img_labeled.height
    imgs_labeled.append(img_labeled)

# Combine vertically
combined_img = Image.new("RGBA", (imgs_labeled[0].width, total_height), (255, 255, 255, 255))
height_added = 0
for img_labeled in imgs_labeled:
    combined_img.paste(img_labeled, (0, height_added))
    height_added += img_labeled.height

# Save final image
combined_img.save(
    f"{RESULTS_DIR}/figure4.png", 
    dpi=(600, 600),
)


# =========================================================================
# SUPP FIGURE 2 -- PLOT AGE DISTRIBUTIONS
# =========================================================================
fig, axs = plt.subplots(
    nrows=len(COHORTS),
    figsize=(10, 6 * len(COHORTS)),
)

if len(COHORTS) == 1:
    axs = [axs]

for i, cohort in enumerate(COHORTS):
    # Read and reformat metadata
    df = pd.read_csv(
        f"{PREPROCESSING_DIR}/{cohort}/{cohort}_promoters_covariates_merged.txt", 
        sep="\t",
    )[["sample", "age", "sex_binary"]]
    df["sex"] = np.where(df["sex_binary"] == 1, "Male", "Female")
    df.drop(columns=["sex_binary"], inplace=True)

    # Bin age in 10-year increments
    bins = list(range(0, 101, 10))
    age_bin_labels = [f"{i}–{i+9}" for i in range(0, 100, 10)]

    df["AGE_BIN"] = pd.cut(
        df["age"],
        bins=bins,
        labels=age_bin_labels,
        include_lowest=True,
        right=False           # intervals are left-inclusive, right-exclusive [0, 10)
    )
    age_counts = df.groupby(["AGE_BIN", "sex"], observed=False)["sample"].count()
    age_counts_pivot = age_counts.unstack("sex").fillna(0)

    # Plot age bins
    ax = axs[i]
    age_counts_pivot.plot(kind="bar", ax=ax, color=["#B359BF", "#4D26BF"])
    for container in ax.containers:
        ax.bar_label(container, label_type="edge", fontsize=8)
    ax.set_ylabel("Count")
    ax.set_title(f"{cohort} Age Distribution by Sex")
    ax.legend(title="Sex")

    ax.set_xlabel("Age range (years)")
    ax.set_xticks(np.arange(len(age_bin_labels)))
    ax.set_xticklabels(age_bin_labels, rotation=45, ha="right")

plt.tight_layout()
plt.savefig(
    f"{RESULTS_DIR}/supp_fig2.png", 
    dpi=600,
)


# =========================================================================
# SUPP TABLE 1 -- BETA VALUES FOR EACH MODEL
# SUPP TABLE 2 -- FEATURE OVERLAPS
# =========================================================================
num_overlap = {}
for cohort in COHORTS:
    # Initialize list for current cohort
    num_overlap[cohort] = []

    # Find coefficients and corresponding promoters from the mdoel
    model = joblib.load(f"{RESULTS_DIR}/{cohort}_promoter/model.joblib")
    coefficients = np.insert(model.coef_, 0, model.intercept_)
    features = pd.read_csv(
        f"{PREPROCESSING_DIR}/{cohort}/{cohort}_train_promoter.tsv", 
        sep="\t",
    ).columns.values[1:]
    features = np.insert(features, 0, "Intercept")

    df_model_betas = pd.DataFrame({
        "Feature": features,
        "Beta": coefficients,
    })

    # Save only features with nonzero coefficients
    df_model_betas = df_model_betas[df_model_betas.Beta != 0]
    df_model_betas.to_csv(
        f"{RESULTS_DIR}/supp_tab1_{cohort}.tsv",
        sep="\t",
        index=False,
    )

    # Find total number of promoters per cohort and overlap with previous clock models
    cohort_promoters = [feature for feature in df_model_betas.Feature.values if not feature in ["Intercept", "sex_binary", "PMI"]]
    num_overlap[cohort].append(len(cohort_promoters))
    for clock_name in ["clock1", "clock2", "clock3"]:
        clock_promoters = pd.read_csv(
            f"{PREPROCESSING_DIR}/{cohort}/{cohort}_{clock_name.upper()}_PROMOTER_SUBSET.txt",
            sep="\t",
        ).columns.values[1:]
        num_overlap[cohort].append(len([promoter for promoter in clock_promoters if promoter in cohort_promoters]))

df_overlap = pd.DataFrame(num_overlap, index=["Total","Lu Clock 1","Lu Clock 2","Lu Clock 2"])
df_overlap.to_csv(
        f"{RESULTS_DIR}/supp_tab2.tsv",
    sep="\t",
)
