## =========================================================
##  Bulk RNA-seq Downstream Analysis
##  PCA and Pathway Enrichment
##  Figure 8p,r and Figure 9d-f
##
##  DEG analysis was performed using CLC Genomics Workbench.
##  This script loads CLC output for downstream analysis.
##
##  Input:  CLC DEG tables, expression matrices
##  Output: PCA plots, pathway enrichment results
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(enrichR)
})

## ---------- Load CLC Output ----------
# iPSC-derived microglia (GEO: GSE325368)
deg_iMG_gd  <- read.csv("data/bulk_rnaseq/iMG_DEG_GD_vs_Control.csv", row.names = 1)
deg_iMG_ma5 <- read.csv("data/bulk_rnaseq/iMG_DEG_MA5_vs_GD.csv", row.names = 1)
iMG_expr    <- read.csv("data/bulk_rnaseq/iMG_bulk_expression.csv", row.names = 1)

# Mouse liver (GEO: GSE325367)
deg_liver_gd  <- read.csv("data/bulk_rnaseq/liver_DEG_GD_vs_Control.csv", row.names = 1)
deg_liver_ma5 <- read.csv("data/bulk_rnaseq/liver_DEG_MA5_vs_GD.csv", row.names = 1)
liver_tpm     <- read.csv("data/bulk_rnaseq/liver_tpm.csv", row.names = 1)

## ---------- PCA (Figure 8p, 9d) ----------
group_colors <- c("Control" = "#95C9FF", "Gaucher" = "#FF9B8E", "Gaucher+MA5" = "#A3E6A4")

create_pca_plot <- function(expr_matrix, sample_groups, title) {
  pca_res <- prcomp(t(log2(expr_matrix + 1)), scale. = TRUE)
  pca_df  <- data.frame(pca_res$x[, 1:2], Group = sample_groups)
  var_exp <- summary(pca_res)$importance[2, 1:2] * 100

  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 4, stroke = 1) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    scale_color_manual(values = group_colors) +
    labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         title = title) +
    theme_bw(base_family = "Arial") +
    theme(legend.position = "right")
}

## ---------- Rescue Gene Identification ----------
deg_params <- read.delim("config/deg_filter_parameters.txt", header = TRUE)

identify_rescue_genes <- function(deg_gd, deg_ma5, params) {
  up_in_gd    <- rownames(deg_gd[deg_gd$logFC > params$logfc_thresh &
                                  deg_gd$adj.P.Val < params$p_thresh, ])
  down_by_ma5 <- rownames(deg_ma5[deg_ma5$logFC < -params$logfc_thresh &
                                   deg_ma5$adj.P.Val < params$p_thresh, ])
  intersect(up_in_gd, down_by_ma5)
}

rescue_genes_iMG   <- identify_rescue_genes(deg_iMG_gd,   deg_iMG_ma5,   deg_params)
rescue_genes_liver <- identify_rescue_genes(deg_liver_gd, deg_liver_ma5, deg_params)

## ---------- Enrichr (Figure 8r, 9e-f) ----------
enrichr_databases <- c("MSigDB_Hallmark_2020", "Reactome_2022",
                        "TRRUST_Transcription_Factors_2019",
                        "BioPlanet_2019", "GO_Biological_Process_2023")

enrichr_iMG   <- enrichr(rescue_genes_iMG,   databases = enrichr_databases)
enrichr_liver <- enrichr(rescue_genes_liver, databases = enrichr_databases)
