## =========================================================
##  snRNA-seq Preprocessing and Cell Type Annotation
##  Figure 3b-d
##
##  Input:  CellRanger output (10x Genomics)
##  Output: combined_seurat_object.rds
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

## ---------- Load CellRanger Output ----------
ctrl.data    <- Read10X(data.dir = "data/Ctrl/outs/filtered_feature_bc_matrix")
gaucher.data <- Read10X(data.dir = "data/Gaucher/outs/filtered_feature_bc_matrix")
ma5.data     <- Read10X(data.dir = "data/Gaucher_MA5/outs/filtered_feature_bc_matrix")

ctrl    <- CreateSeuratObject(counts = ctrl.data,    project = "Control")
gaucher <- CreateSeuratObject(counts = gaucher.data, project = "Gaucher")
ma5     <- CreateSeuratObject(counts = ma5.data,     project = "Gaucher+MA5")

## ---------- Quality Control ----------
qc_params <- read.delim("config/qc_parameters.txt", header = TRUE)

ctrl[["percent.mt"]]    <- PercentageFeatureSet(ctrl,    pattern = "^mt-")
gaucher[["percent.mt"]] <- PercentageFeatureSet(gaucher, pattern = "^mt-")
ma5[["percent.mt"]]     <- PercentageFeatureSet(ma5,     pattern = "^mt-")

filter_cells <- function(obj, params) {
  subset(obj,
    subset = percent.mt < params$mt_max &
             nFeature_RNA > params$nFeature_min &
             nFeature_RNA < params$nFeature_max)
}
ctrl    <- filter_cells(ctrl,    qc_params)
gaucher <- filter_cells(gaucher, qc_params)
ma5     <- filter_cells(ma5,     qc_params)

## ---------- Integration and Clustering ----------
cluster_params <- read.delim("config/clustering_parameters.txt", header = TRUE)

combined <- merge(ctrl, y = c(gaucher, ma5))
combined <- NormalizeData(combined) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

combined <- FindNeighbors(combined, dims = 1:cluster_params$n_dims)
combined <- FindClusters(combined, resolution = cluster_params$resolution)
combined <- RunUMAP(combined, dims = 1:cluster_params$n_dims)

## ---------- Cell Type Annotation ----------
annotation_markers <- read.delim("config/celltype_markers.txt", header = TRUE)

celltype_colors <- c(
  "Neuron" = "#3498DB", "Astrocyte" = "#E67E22", "Oligodendrocyte" = "#9B59B6",
  "Microglia" = "#E74C3C", "OPC" = "#F1C40F", "Endothelial cell" = "#1ABC9C"
)

## ---------- Marker Heatmap (Figure 3b) ----------
all_markers <- FindAllMarkers(combined, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.5)
top5 <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(combined, features = top5$gene) + NoLegend()

## ---------- UMAP Visualization (Figure 3c-d) ----------
DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined, reduction = "umap", group.by = "cell_type",
        cols = celltype_colors, label = FALSE)

saveRDS(combined, file = "output/combined_seurat_object.rds")
