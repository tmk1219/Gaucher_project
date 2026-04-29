#!/usr/bin/env python3
"""
Microglia Subcluster UMAP with IFN/DAM Score Overlay
Figure 5a-c

Input:  microglia_data_with_all_scores.csv
Output: 9-panel UMAP (3 conditions x [cluster, IFN, DAM])
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
data = pd.read_csv("data/microglia/microglia_data_with_all_scores.csv")

## ---------- Group Mapping ----------
group_map = {"Control": "Control", "Gaucher": "GD", "Gaucher+MA5": "GD+MA-5"}
data["Group"] = data["orig.ident"].map(group_map)
group_order = ["Control", "GD", "GD+MA-5"]

## ---------- Cluster Colors ----------
cluster_colors = {
    0: "#E74C3C", 1: "#3498DB", 2: "#2ECC71", 3: "#F39C12", 4: "#9B59B6"
}

## ---------- Score Colormaps ----------
score_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "score", ["#F0F0F0", "#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"])

## ---------- 9-Panel Figure ----------
import os
os.makedirs("output/microglia_umaps", exist_ok=True)

fig, axes = plt.subplots(3, 3, figsize=(15, 14))
panel_labels = ["Cluster", "IFN score", "DAM score"]
score_cols = [None, "Score_IFN", "Score_DAM"]

for col_idx, (label, score_col) in enumerate(zip(panel_labels, score_cols)):
    for row_idx, group in enumerate(group_order):
        ax = axes[row_idx, col_idx]
        subset = data[data["Group"] == group]

        if score_col is None:
            # Cluster UMAP
            for cl in sorted(subset["seurat_clusters"].unique()):
                cl_data = subset[subset["seurat_clusters"] == cl]
                ax.scatter(cl_data["UMAP_1"], cl_data["UMAP_2"],
                          c=[cluster_colors.get(cl, "grey")],
                          s=3, alpha=0.7, label=str(cl), rasterized=True)
            if row_idx == 0:
                ax.legend(title="Cluster", markerscale=3, fontsize=8,
                         loc="upper right", framealpha=0.8)
        else:
            # Score UMAP
            vmin = data[score_col].quantile(0.01)
            vmax = data[score_col].quantile(0.99)
            ax.scatter(subset["UMAP_1"], subset["UMAP_2"],
                      c=subset[score_col], cmap=score_cmap,
                      s=3, alpha=0.7, vmin=vmin, vmax=vmax, rasterized=True)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")

        if col_idx == 0:
            ax.set_ylabel(group, fontsize=14, fontweight="bold")
        if row_idx == 0:
            ax.set_title(label, fontsize=14, fontweight="bold")

plt.tight_layout()
fig.savefig("output/microglia_umaps/Fig5_9panel.pdf", bbox_inches="tight", dpi=300)
fig.savefig("output/microglia_umaps/Fig5_9panel.png", bbox_inches="tight", dpi=300)
plt.close()
