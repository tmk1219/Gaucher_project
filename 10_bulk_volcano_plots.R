## =========================================================
##  Volcano Plot Visualization
##  Supplementary
##
##  Input:  CLC DEG tables
##  Output: Volcano plots (PDF/PNG)
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

## ---------- Load DEG Tables ----------
deg_iMG_gd    <- read.csv("data/bulk_rnaseq/iMG_DEG_GD_vs_Control.csv", row.names = 1)
deg_iMG_ma5   <- read.csv("data/bulk_rnaseq/iMG_DEG_MA5_vs_GD.csv", row.names = 1)
deg_liver_gd  <- read.csv("data/bulk_rnaseq/liver_DEG_GD_vs_Control.csv", row.names = 1)
deg_liver_ma5 <- read.csv("data/bulk_rnaseq/liver_DEG_MA5_vs_GD.csv", row.names = 1)

## ---------- Volcano Plot Function ----------
colors <- c("Up" = "#FF7A6E", "Down" = "#5BA3E0", "NS" = "grey80")

create_volcano <- function(deg_table, title, fc_thresh = 1, p_thresh = 0.05) {
  deg_table$gene <- rownames(deg_table)
  deg_table$direction <- ifelse(
    deg_table$logFC > fc_thresh & deg_table$adj.P.Val < p_thresh, "Up",
    ifelse(deg_table$logFC < -fc_thresh & deg_table$adj.P.Val < p_thresh, "Down", "NS"))

  top_genes <- deg_table %>%
    filter(direction != "NS") %>%
    arrange(adj.P.Val) %>%
    head(20)

  ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    geom_point(size = 1, alpha = 0.7) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3, max.overlaps = 15) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(adj. P)") +
    theme_bw(base_family = "Arial")
}

## ---------- Generate Plots ----------
dir.create("output/volcano", showWarnings = FALSE, recursive = TRUE)

ggsave("output/volcano/iMG_GD_vs_Ctrl.pdf",
       create_volcano(deg_iMG_gd, "iPSC Microglia: GD vs Control"), width = 8, height = 6)
ggsave("output/volcano/iMG_MA5_vs_GD.pdf",
       create_volcano(deg_iMG_ma5, "iPSC Microglia: GD+MA-5 vs GD"), width = 8, height = 6)
ggsave("output/volcano/liver_GD_vs_Ctrl.pdf",
       create_volcano(deg_liver_gd, "Liver: GD vs Control"), width = 8, height = 6)
ggsave("output/volcano/liver_MA5_vs_GD.pdf",
       create_volcano(deg_liver_ma5, "Liver: GD+MA-5 vs GD"), width = 8, height = 6)
