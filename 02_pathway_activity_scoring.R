## =========================================================
##  Pathway Activity Scoring (12 MSigDB Hallmark Pathways)
##  Figure 3e
##
##  Input:  combined_seurat_object.rds
##  Output: pathway delta score heatmap
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

combined <- readRDS("output/combined_seurat_object.rds")

## ---------- Load Pathway Gene Sets ----------
pathway_gene_sets <- readRDS("config/pathway_gene_sets.rds")

## ---------- Compute Per-Cell Pathway Scores ----------
for (pw_name in names(pathway_gene_sets)) {
  genes <- pathway_gene_sets[[pw_name]]
  genes_present <- intersect(genes, rownames(combined))
  message(sprintf("  %s: %d / %d genes", pw_name, length(genes_present), length(genes)))
  combined <- AddModuleScore(combined, features = list(genes_present),
                              name = paste0("Score_", pw_name))
}

## ---------- Delta Score Calculation ----------
# Delta = (condition mean - Control mean) x 10
compute_delta_scores <- function(seurat_obj, score_cols, cell_types) {
  meta <- seurat_obj@meta.data
  ctrl_cells <- meta$orig.ident == "Control"

  results <- list()
  for (ct in cell_types) {
    ct_cells <- meta$cell_type == ct
    for (sc in score_cols) {
      ctrl_mean <- mean(meta[ctrl_cells & ct_cells, sc], na.rm = TRUE)
      for (grp in unique(meta$orig.ident)) {
        grp_cells <- meta$orig.ident == grp & ct_cells
        delta <- (mean(meta[grp_cells, sc], na.rm = TRUE) - ctrl_mean) * 10
        results[[length(results) + 1]] <- data.frame(
          CellType = ct, Pathway = sc, Group = grp, Delta = delta)
      }
    }
  }
  do.call(rbind, results)
}

cell_types <- c("Microglia", "Neuron", "Astrocyte",
                "Oligodendrocyte", "OPC", "Endothelial cell")
score_cols <- grep("^Score_", colnames(combined@meta.data), value = TRUE)
delta_df <- compute_delta_scores(combined, score_cols, cell_types)

saveRDS(combined, file = "output/combined_seurat_with_scores.rds")
write.csv(delta_df, "output/pathway_delta_scores.csv", row.names = FALSE)
