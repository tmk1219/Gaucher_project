## =========================================================
##  Spearman Correlation Heatmap
##  Figure 7b (correlation panel)
##
##  Input:  Pre-computed correlation matrices
##  Output: Correlation heatmap (PDF/PNG)
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
})

## ---------- Load Data ----------
cor_matrix  <- read.csv("data/gsea/Correlation_table_genes_Spearman.csv",
                         row.names = 1, check.names = FALSE)
padj_matrix <- read.csv("data/gsea/Correlation_table_genes_Spearman_padjusted.csv",
                         row.names = 1, check.names = FALSE)

## ---------- Significance Labels ----------
get_sig_label <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",
  ifelse(p < 0.1,   "^", ""))))
}

## ---------- Heatmap ----------
melted <- melt(as.matrix(cor_matrix))
colnames(melted) <- c("Row", "Col", "Value")
melted$Label <- get_sig_label(as.vector(as.matrix(padj_matrix)))

ggplot(melted, aes(Col, Row, fill = Value)) +
  geom_tile(colour = "black", linewidth = 0.1) +
  scale_fill_gradient2(high = "red3", low = "blue3", mid = "white",
                        limits = c(-1, 1), oob = squish) +
  geom_text(aes(label = Label), size = 4, color = "white", fontface = "bold") +
  labs(x = "", y = "", fill = expression(Spearman ~ rho)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
