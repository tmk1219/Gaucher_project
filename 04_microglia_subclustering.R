## =========================================================
##  Microglia Subclustering Analysis
##  Figure 5a-c
##
##  Input:  combined_seurat_all_signatures.rds
##  Output: microglia_subclustered.rds
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

combined <- readRDS("output/combined_seurat_all_signatures.rds")

## ---------- Parameters ----------
mg_params <- read.delim("config/microglia_clustering_parameters.txt", header = TRUE)

## ---------- Extract Microglia ----------
Idents(combined) <- "cell_type"
microglia <- subset(combined, idents = "Microglia")

# Remove MRC1+ macrophage contamination
microglia <- subset(microglia, subset = Mrc1 < mg_params$mrc1_threshold)

## ---------- Re-clustering ----------
microglia <- NormalizeData(microglia) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

microglia <- FindNeighbors(microglia, dims = 1:mg_params$n_dims)
microglia <- FindClusters(microglia, resolution = mg_params$resolution)
microglia <- RunUMAP(microglia, dims = 1:mg_params$n_dims)

## ---------- Cluster Proportions ----------
prop_table <- table(microglia$seurat_clusters, microglia$orig.ident) %>%
  prop.table(margin = 2) %>% round(3)
print(prop_table)
write.csv(as.data.frame.matrix(prop_table), "output/microglia_cluster_proportions.csv")

## ---------- UMAP (Figure 5a) ----------
DimPlot(microglia, reduction = "umap", split.by = "orig.ident",
        label = TRUE, repel = TRUE)

## ---------- Score Visualization (Figure 5b-c) ----------
FeaturePlot(microglia, features = "Score_IFN1",
            split.by = "orig.ident", cols = c("grey90", "red3"))
FeaturePlot(microglia, features = "Score_DAM1",
            split.by = "orig.ident", cols = c("grey90", "red3"))

saveRDS(microglia, file = "output/microglia_subclustered.rds")
