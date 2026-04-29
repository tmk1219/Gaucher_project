#!/usr/bin/env python3
"""
Signature Score UMAP Visualization (Python)
Figure 4a-e (UMAP panels)

Generates 3-condition UMAP panels for each signature score:
  DAM, IFN, Microglia Activation, NDMG, SASP

Input:  all_cells_scores_data.csv
Output: Score UMAP PNG/PDF per signature
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams.update({
    'font.family': 'Arial',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

## ---------- Load Data ----------
data = pd.read_csv("data/scores/all_cells_scores_data.csv")

## ---------- Score Columns ----------
score_info = {
    "Score_DAM":                  "DAM score",
    "Score_IFN":                  "IFN score",
    "Score_MicrogliaActivation":  "Microglia activation score",
    "Score_NDMG":                 "NDMG score",
    "Score_SASP":                 "SASP score",
}

## ---------- Group Mapping ----------
group_map = {"Control": "Control", "Gaucher": "GD", "Gaucher+MA5": "GD+MA-5"}
data["Group"] = data["orig.ident"].map(group_map)
group_order = ["Control", "GD", "GD+MA-5"]

## ---------- Color Palette ----------
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "score_cmap", ["#F0F0F0", "#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"])

## ---------- Generate Plots ----------
import os
os.makedirs("output/score_umaps", exist_ok=True)

for score_col, score_label in score_info.items():
    if score_col not in data.columns:
        print(f"  Skipping {score_col}: not found")
        continue

    vmin = data[score_col].quantile(0.01)
    vmax = data[score_col].quantile(0.99)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for ax, group in zip(axes, group_order):
        subset = data[data["Group"] == group]
        sc = ax.scatter(
            subset["UMAP_1"], subset["UMAP_2"],
            c=subset[score_col], cmap=cmap,
            s=1, alpha=0.8, vmin=vmin, vmax=vmax,
            rasterized=True
        )
        ax.set_title(group, fontsize=14, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")

        # Small axis arrows
        x0, y0 = ax.get_xlim()[0], ax.get_ylim()[0]
        length = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.15
        ax.annotate("", xy=(x0 + length, y0), xytext=(x0, y0),
                     arrowprops=dict(arrowstyle="->", lw=1.5))
        ax.annotate("", xy=(x0, y0 + length), xytext=(x0, y0),
                     arrowprops=dict(arrowstyle="->", lw=1.5))

    fig.suptitle(score_label, fontsize=16, fontweight="bold", y=1.02)
    plt.tight_layout()

    fig.savefig(f"output/score_umaps/{score_col}_3groups.pdf",
                bbox_inches="tight", dpi=300)
    fig.savefig(f"output/score_umaps/{score_col}_3groups.png",
                bbox_inches="tight", dpi=300)
    plt.close()
    print(f"  {score_label}: done")
