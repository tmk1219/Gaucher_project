#!/usr/bin/env python3
"""
SCAFE Enhancer/Promoter Ranking and Visualization
Figure 6b,f (ranking tables), 6d,e (arch diagrams)

Input:  SCAFE CRE activity data, Cicero co-accessibility
Output: Ranking tables, enhancer-promoter arch plots (PDF/PNG)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams.update({
    'font.family': 'Arial',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

os.makedirs("output/scafe_figures", exist_ok=True)

## ---------- Load CRE Activity Data ----------
promoter_activity = pd.read_csv("data/scafe/promoter_activity.csv", index_col=0)
enhancer_activity = pd.read_csv("data/scafe/enhancer_activity.csv", index_col=0)

## ---------- Rescue Score Calculation ----------
def compute_rescue_score(activity_df, pattern="GD_up_MA5_down"):
    """Rank promoters/enhancers by GD activation and MA-5 rescue."""
    ctrl = activity_df["Control"]
    gd   = activity_df["GD"]
    ma5  = activity_df["MA5"]

    if pattern == "GD_up_MA5_down":
        fc = gd - ctrl
        rescue_pct = (gd - ma5) / (gd - ctrl + 1e-10) * 100
    else:  # GD_down_MA5_up
        fc = ctrl - gd
        rescue_pct = (ma5 - gd) / (ctrl - gd + 1e-10) * 100

    # Rescue weight
    weight = np.where(rescue_pct >= 80, 1.5,
             np.where(rescue_pct >= 60, 1.0,
             np.where(rescue_pct >= 40, 0.5, 0.2)))

    score = fc * weight
    activity_df["rescue_score"] = score
    activity_df["rescue_pct"]   = rescue_pct
    return activity_df.sort_values("rescue_score", ascending=False)

## ---------- Top 10 Rankings (Figure 6b, 6f) ----------
promoter_ranked = compute_rescue_score(promoter_activity, "GD_up_MA5_down")
enhancer_ranked = compute_rescue_score(enhancer_activity, "GD_down_MA5_up")

promoter_ranked.head(10).to_csv("output/scafe_figures/top10_promoters_GD_up.csv")
enhancer_ranked.head(10).to_csv("output/scafe_figures/top10_enhancers_GD_down.csv")

## ---------- Arch Diagram Function ----------
def plot_enhancer_promoter_arch(gene, promoter_pos, enhancer_pos, co_activity,
                                 conditions=["Control", "GD", "GD+MA-5"]):
    """Draw enhancer-promoter connection arch with co-activity scores."""
    fig, ax = plt.subplots(figsize=(10, 3))

    # Genomic region
    region_start = min(promoter_pos, enhancer_pos) - 5000
    region_end   = max(promoter_pos, enhancer_pos) + 5000

    # Promoter and enhancer positions
    ax.axvline(promoter_pos, color="blue", linewidth=2, label="Promoter")
    ax.axvline(enhancer_pos, color="red",  linewidth=2, label="Enhancer")

    # Arch
    mid = (promoter_pos + enhancer_pos) / 2
    width = abs(enhancer_pos - promoter_pos)
    height = width * 0.3

    for i, (cond, co_act) in enumerate(zip(conditions, co_activity)):
        alpha = min(co_act * 5, 1.0)
        arc = patches.Arc((mid, 0), width, height * 2,
                          angle=0, theta1=0, theta2=180,
                          color=["#6BAED6", "#FC9272", "#74C476"][i],
                          linewidth=2, alpha=alpha)
        ax.add_patch(arc)
        ax.text(mid, height * (1.1 + i * 0.15),
                f"{cond}: {co_act:.3f}", ha="center", fontsize=9)

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.5, height * 2)
    ax.set_title(gene, fontsize=14, fontweight="bold")
    ax.set_xlabel("Genomic position", fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    plt.tight_layout()
    return fig

## ---------- Tbk1 (Figure 6d) ----------
fig = plot_enhancer_promoter_arch(
    "Tbk1",
    promoter_pos=121586913,
    enhancer_pos=121600183,
    co_activity=[0.104, 0.189, 0.041]
)
fig.savefig("output/scafe_figures/Fig6d_Tbk1_arch.pdf", bbox_inches="tight")
plt.close()

## ---------- Ifit Cluster (Figure 6e) ----------
# Shared enhancer at chr19:34,654,945-34,655,446
# Co-activity with all 4 Ifit gene promoters (GD-specific)
ifit_genes = {
    "Ifit1": {"promoter": 34627000, "co_activity": [0, 0.045, 0]},
    "Ifit2": {"promoter": 34637000, "co_activity": [0, 0.019, 0]},
    "Ifit3": {"promoter": 34645000, "co_activity": [0, 0.032, 0]},
    "Ifit3b": {"promoter": 34650000, "co_activity": [0, 0.025, 0]},
}
shared_enhancer = 34655195

for gene, info in ifit_genes.items():
    fig = plot_enhancer_promoter_arch(
        gene, info["promoter"], shared_enhancer, info["co_activity"])
    fig.savefig(f"output/scafe_figures/Fig6e_{gene}_arch.pdf", bbox_inches="tight")
    plt.close()
