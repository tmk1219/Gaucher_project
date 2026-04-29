## =========================================================
##  Microglia Cluster Marker Genes and Enrichment Analysis
##  Figure 5d-g
##
##  Input:  microglia_subclustered.rds
##  Output: cluster marker CSVs, Enrichr results
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(enrichR)
})

microglia <- readRDS("output/microglia_subclustered.rds")
Idents(microglia) <- "seurat_clusters"

## ---------- Find Cluster Markers ----------
all_markers <- FindAllMarkers(microglia, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.5,
                               test.use = "wilcox")

## ---------- Export Per-Cluster Markers ----------
for (cluster_id in unique(all_markers$cluster)) {
  markers <- all_markers %>%
    filter(cluster == cluster_id) %>%
    arrange(desc(avg_log2FC))
  write.csv(markers, sprintf("output/cluster_%s_markers.csv", cluster_id),
            row.names = FALSE)
}

## ---------- Enrichr Pathway Analysis ----------
enrichr_databases <- c("GO_Biological_Process_2023", "KEGG_2021_Human")

for (cluster_id in unique(all_markers$cluster)) {
  top_genes <- all_markers %>%
    filter(cluster == cluster_id) %>%
    arrange(desc(avg_log2FC)) %>%
    head(250) %>%
    pull(gene)

  if (length(top_genes) > 0) {
    results <- enrichr(top_genes, databases = enrichr_databases)
    for (db in names(results)) {
      write.csv(results[[db]],
                sprintf("output/enrichr_cluster%s_%s.csv", cluster_id, db),
                row.names = FALSE)
    }
  }
}

## ---------- Combined Clusters 3+4 (Figure 5g) ----------
cluster34_markers <- all_markers %>%
  filter(cluster %in% c("3", "4")) %>%
  arrange(desc(avg_log2FC))
write.csv(cluster34_markers, "output/cluster_3_4_combined_markers.csv", row.names = FALSE)

enrichr_34 <- enrichr(head(cluster34_markers$gene, 250), databases = enrichr_databases)
for (db in names(enrichr_34)) {
  write.csv(enrichr_34[[db]],
            sprintf("output/enrichr_cluster34_%s.csv", db), row.names = FALSE)
}
