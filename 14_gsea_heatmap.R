## =========================================================
##  GSEA NES Heatmap
##  Figure 7b (GSEA panel), 8q, 9g
##
##  Input:  Pre-computed GSEA NES and p-adjusted matrices
##  Output: GSEA heatmap (PDF/PNG)
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
})

## ---------- Load Data ----------
nes_matrix <- read.csv("data/gsea/GSEA_matrix_NES.csv",
                        row.names = 1, check.names = FALSE)
nes_padj   <- read.csv("data/gsea/GSEA_matrix_p.adjusted.csv",
                        row.names = 1, check.names = FALSE)

## ---------- Significance Labels ----------
get_sig_label <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",
  ifelse(p < 0.1,   "^", ""))))
}

label_mat <- apply(as.matrix(nes_padj), c(1, 2), get_sig_label)

## ---------- Heatmap ----------
nes_clamped <- pmin(pmax(as.matrix(nes_matrix), -2), 2)

melted <- melt(nes_clamped)
colnames(melted) <- c("Pathway", "Dataset", "NES")
melted$Label <- as.vector(label_mat)

ggplot(melted, aes(Dataset, Pathway, fill = NES)) +
  geom_tile(colour = "black", linewidth = 0.1) +
  scale_fill_gradient2(high = "red3", low = "blue3", mid = "white",
                        limits = c(-2, 2), oob = squish) +
  geom_text(aes(label = Label), size = 4, color = "white", fontface = "bold") +
  labs(x = "", y = "", fill = "NES") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
