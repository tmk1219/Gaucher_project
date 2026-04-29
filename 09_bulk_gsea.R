## =========================================================
##  GSEA Analysis on Bulk RNA-seq DEG Results
##  Figure 8q, 9g
##
##  Input:  CLC DEG tables
##  Output: GSEA NES matrices, pathway heatmaps
## =========================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
})

## ---------- Load DEG Tables ----------
deg_iMG_gd    <- read.csv("data/bulk_rnaseq/iMG_DEG_GD_vs_Control.csv", row.names = 1)
deg_iMG_ma5   <- read.csv("data/bulk_rnaseq/iMG_DEG_MA5_vs_GD.csv", row.names = 1)
deg_liver_gd  <- read.csv("data/bulk_rnaseq/liver_DEG_GD_vs_Control.csv", row.names = 1)
deg_liver_ma5 <- read.csv("data/bulk_rnaseq/liver_DEG_MA5_vs_GD.csv", row.names = 1)

## ---------- Pathway Databases ----------
m_hallmark <- msigdbr(species = "Mus musculus", collection = "H")
m_c2       <- msigdbr(species = "Mus musculus", collection = "C2")
m_kegg     <- m_c2[grepl("^KEGG_", m_c2$gs_name), ]
m_reactome <- m_c2[grepl("^REACTOME_", m_c2$gs_name), ]

pathways <- split(
  c(m_hallmark$gene_symbol, m_kegg$gene_symbol, m_reactome$gene_symbol),
  c(m_hallmark$gs_name, m_kegg$gs_name, m_reactome$gs_name)
)

## ---------- Run GSEA ----------
run_gsea <- function(deg_table) {
  ranked <- setNames(deg_table$logFC, rownames(deg_table))
  ranked <- sort(ranked, decreasing = TRUE)
  fgseaMultilevel(pathways, ranked, nPermSimple = 10000)
}

gsea_iMG_gd    <- run_gsea(deg_iMG_gd)
gsea_iMG_ma5   <- run_gsea(deg_iMG_ma5)
gsea_liver_gd  <- run_gsea(deg_liver_gd)
gsea_liver_ma5 <- run_gsea(deg_liver_ma5)

## ---------- Export Results ----------
dir.create("output/gsea", showWarnings = FALSE)
write.csv(as.data.frame(gsea_iMG_gd),    "output/gsea/gsea_iMG_GD_vs_Ctrl.csv")
write.csv(as.data.frame(gsea_iMG_ma5),   "output/gsea/gsea_iMG_MA5_vs_GD.csv")
write.csv(as.data.frame(gsea_liver_gd),  "output/gsea/gsea_liver_GD_vs_Ctrl.csv")
write.csv(as.data.frame(gsea_liver_ma5), "output/gsea/gsea_liver_MA5_vs_GD.csv")
