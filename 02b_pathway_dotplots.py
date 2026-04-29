#!/usr/bin/env python3
"""
Pathway Enrichment Dot Plots
Figure 3f and Figure 5d-g

Input:  Enrichr output CSVs
Output: Dot plots (PDF/PNG)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams.update({
    'font.family': 'Arial',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

## ---------- Load Enrichr Results ----------
# Figure 3f: Layer V neuron pathway analysis
fig3f_concave = pd.read_csv("data/pathway/Figure_3F/L5_Concave_R.csv")
fig3f_convex  = pd.read_csv("data/pathway/Figure_3F/L5_Convex_R.csv")

# Figure 5d-g: Microglia cluster pathway analysis
fig5_data = {}
for cluster in ["0", "1", "2", "3_4"]:
    for db in ["GO_Biological_Process_2023", "KEGG_2021_Human"]:
        fname = f"output/enrichr_cluster{cluster}_{db}.csv"
        if os.path.exists(fname):
            fig5_data[f"cluster{cluster}_{db}"] = pd.read_csv(fname)

## ---------- Dot Plot Function ----------
def create_pathway_plot(df, title, n_terms=10, figsize=(8, 5),
                        size_col="Odds Ratio", color_col="Adjusted P-value"):
    # Clean column names
    df.columns = [c.replace('\ufeff', '') for c in df.columns]

    # Sort by p-value and take top terms
    df = df.sort_values(color_col).head(n_terms)
    df = df.iloc[::-1]  # reverse for bottom-to-top

    # Clean term names
    df["Term_clean"] = df["Term"].str.split(" \\(GO").str[0]
    df["Term_clean"] = df["Term_clean"].str[:60]

    # -log10 p-value for color
    df["neg_log10_p"] = -np.log10(df[color_col].clip(lower=1e-20))

    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        df["neg_log10_p"], range(len(df)),
        s=df[size_col] * 50,
        c=df["neg_log10_p"],
        cmap="RdBu_r",
        edgecolors="black", linewidths=0.5,
        alpha=0.8
    )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term_clean"], fontsize=10)
    ax.set_xlabel(r"$-\log_{10}$(Adjusted P-value)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    plt.colorbar(scatter, ax=ax, label=r"$-\log_{10}$(P)", shrink=0.6)
    plt.tight_layout()
    return fig

## ---------- Generate Figure 3f ----------
os.makedirs("output/pathway_plots", exist_ok=True)

fig = create_pathway_plot(fig3f_concave, "Layer V Neurons: GD up / MA-5 down",
                          n_terms=10, figsize=(10, 5))
fig.savefig("output/pathway_plots/Fig3f_concave.pdf", bbox_inches="tight")
plt.close()

fig = create_pathway_plot(fig3f_convex, "Layer V Neurons: GD down / MA-5 up (SynGO)",
                          n_terms=10, figsize=(10, 5))
fig.savefig("output/pathway_plots/Fig3f_convex.pdf", bbox_inches="tight")
plt.close()

## ---------- Generate Figure 5d-g ----------
for key, df in fig5_data.items():
    fig = create_pathway_plot(df, key.replace("_", " "), n_terms=10)
    fig.savefig(f"output/pathway_plots/Fig5_{key}.pdf", bbox_inches="tight")
    plt.close()
