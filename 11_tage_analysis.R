## =========================================================
##  Transcriptomic Age (tAge) Analysis
##  Figure 7a,c
##
##  tAge was computed using the tAge framework:
##  https://github.com/Gladyshev-Lab/tAge
##  Tyshkovskiy et al., 2024
##
##  Rodent multi-tissue chronological and mortality clocks
##  (Bayesian Ridge algorithm)
##
##  Input:  combined_seurat_object.rds, clock model
##  Output: tAge predictions per cell type and condition
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(metafor)
})

combined <- readRDS("output/combined_seurat_object.rds")

## ---------- Metacell Construction ----------
metacell_params <- read.delim("config/metacell_parameters.txt", header = TRUE)

construct_metacells <- function(seurat_obj, cell_type, params) {
  cells <- WhichCells(seurat_obj, expression = cell_type_annotation == cell_type)
  counts <- GetAssayData(seurat_obj, slot = "counts")[, cells]

  total_reads <- colSums(counts)
  metacells <- list()
  remaining <- seq_along(cells)

  while (length(remaining) > 0) {
    cumsum_reads <- cumsum(total_reads[remaining])
    idx <- which(cumsum_reads >= params$target_reads)[1]
    if (is.na(idx)) break
    metacells[[length(metacells) + 1]] <- remaining[1:idx]
    remaining <- remaining[(idx + 1):length(remaining)]
  }

  mc_counts <- sapply(metacells, function(idx) {
    rowSums(counts[, idx, drop = FALSE])
  })

  keep <- colSums(mc_counts) >= params$min_reads
  mc_counts[, keep]
}

## ---------- Normalization ----------
norm_params <- read.delim("config/normalization_parameters.txt", header = TRUE)

normalize_metacells <- function(mc_counts, params) {
  pct_expressed <- rowMeans(mc_counts >= params$min_reads_per_gene)
  mc_counts <- mc_counts[pct_expressed >= params$min_pct_expressed, ]

  log_expr <- log2(mc_counts + 1)

  median_expr <- apply(log_expr, 1, median)
  ratios <- sweep(log_expr, 1, median_expr, "-")
  size_factors <- apply(ratios, 2, median)
  normalized <- sweep(log_expr, 2, size_factors, "-")

  t(scale(t(normalized)))
}

## ---------- Apply tAge Clock ----------
# Clock model from: https://github.com/Gladyshev-Lab/tAge
# See Tyshkovskiy et al., 2024 for model training details
predict_tage <- function(scaled_expr, clock_model) {
  clock_genes <- names(clock_model$coefficients)
  common_genes <- intersect(rownames(scaled_expr), clock_genes)

  missing_genes <- setdiff(clock_genes, common_genes)
  if (length(missing_genes) > 0) {
    imputed <- matrix(clock_model$imputation_values[missing_genes],
                      nrow = length(missing_genes), ncol = ncol(scaled_expr))
    rownames(imputed) <- missing_genes
    scaled_expr <- rbind(scaled_expr[common_genes, ], imputed)
  }

  expr_subset <- scaled_expr[clock_genes, ]
  tage <- clock_model$intercept + colSums(expr_subset * clock_model$coefficients[clock_genes])
  return(tage)
}

## ---------- Group Comparison ----------
compare_tage <- function(tage_values, groups) {
  df <- data.frame(tage = tage_values, group = groups)
  pairwise_results <- pairwise.t.test(df$tage, df$group, p.adjust.method = "BH")
  return(pairwise_results)
}

## ---------- Run for All Cell Types ----------
cell_types <- c("Microglia", "Neuron", "Astrocyte",
                "Oligodendrocyte", "OPC", "Endothelial cell")

tage_results <- list()
for (ct in cell_types) {
  message(sprintf("Processing: %s", ct))
  mc <- construct_metacells(combined, ct, metacell_params)
  if (ncol(mc) < 3) { message("  Skipped: too few metacells"); next }
  scaled <- normalize_metacells(mc, norm_params)
  # tage <- predict_tage(scaled, clock_model)  # Requires clock model from Gladyshev-Lab/tAge
  # tage_results[[ct]] <- tage
}
