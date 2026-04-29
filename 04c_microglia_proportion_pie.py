#!/usr/bin/env python3
"""
Microglia Cluster Proportion Pie Charts
Figure 5a (proportion panels)

Input:  microglia_data_with_all_scores.csv
Output: 3 pie charts (Control, GD, GD+MA-5)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

mpl.rcParams.update({
    'font.family': 'Arial',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

## ---------- Load Data ----------
data = pd.read_csv("data/microglia/microglia_data_with_all_scores.csv")

group_map = {"Control": "Control", "Gaucher": "GD", "Gaucher+MA5": "GD+MA-5"}
data["Group"] = data["orig.ident"].map(group_map)
group_order = ["Control", "GD", "GD+MA-5"]

## ---------- Cluster Colors ----------
cluster_colors = ["#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6"]

## ---------- Pie Charts ----------
os.makedirs("output/microglia_umaps", exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

for ax, group in zip(axes, group_order):
    subset = data[data["Group"] == group]
    counts = subset["seurat_clusters"].value_counts().sort_index()
    pcts = counts / counts.sum() * 100

    wedges, texts, autotexts = ax.pie(
        pcts, labels=[f"C{i}" for i in pcts.index],
        colors=[cluster_colors[i] for i in pcts.index],
        autopct=lambda p: f"{p:.1f}%" if p > 3 else "",
        startangle=90, textprops={"fontsize": 10}
    )
    ax.set_title(group, fontsize=14, fontweight="bold")

plt.tight_layout()
fig.savefig("output/microglia_umaps/Fig5_pie_charts.pdf", bbox_inches="tight")
fig.savefig("output/microglia_umaps/Fig5_pie_charts.png", bbox_inches="tight", dpi=300)
plt.close()
