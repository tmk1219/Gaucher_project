## =========================================================
##  Disease Signature Scoring
##  DAM, IFN, Microglia Activation, NDMG, SASP
##  Figure 4a-e
##
##  Input:  combined_seurat_with_scores.rds
##  Output: UMAP and violin plots per signature
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

combined <- readRDS("output/combined_seurat_with_scores.rds")

## ---------- Load Signature Gene Lists ----------
# DAM (278 genes) - Gulen et al., 2023
# IFN (28 genes) - Gulen et al., 2023
# Microglia Activation (7 genes) - Shimizu et al., 2023
# NDMG (116 genes) - Gulen et al., 2023
# SASP (10 genes) - Victorelli et al., 2023
signature_lists <- readRDS("config/signature_gene_lists.rds")

## ---------- Calculate Scores ----------
for (sig_name in names(signature_lists)) {
  genes <- signature_lists[[sig_name]]
  genes_present <- intersect(genes, rownames(combined))
  message(sprintf("  %s: %d / %d genes found", sig_name,
                  length(genes_present), length(genes)))
  combined <- AddModuleScore(combined, features = list(genes_present),
                              name = paste0("Score_", sig_name))
}

## ---------- UMAP Visualization (Figure 4, left panels) ----------
score_cols <- paste0("Score_", names(signature_lists), "1")
for (score in score_cols) {
  if (score %in% colnames(combined@meta.data)) {
    print(FeaturePlot(combined, features = score, split.by = "orig.ident",
                       cols = c("grey90", "red3")))
  }
}

## ---------- Violin Plots by Cell Type (Figure 4, right panels) ----------
group_colors <- c("Control" = "#6BAED6", "GD" = "#FC9272", "GD+MA-5" = "#74C476")
celltype_order <- c("Microglia", "Neuron", "Astrocyte",
                    "Oligodendrocyte", "OPC", "Endothelial cell")

for (score_col in score_cols) {
  if (score_col %in% colnames(combined@meta.data)) {
    VlnPlot(combined, features = score_col, group.by = "orig.ident",
            split.by = "cell_type", cols = group_colors)
  }
}

saveRDS(combined, file = "output/combined_seurat_all_signatures.rds")
