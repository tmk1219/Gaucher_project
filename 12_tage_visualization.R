## =========================================================
##  tAge Visualization: Dotplots and Gene Contribution
##  Figure 7a,c,d
##
##  Input:  Pre-computed tAge predictions, gene contribution tables
##  Output: Dotplots, barplots (PDF/PNG)
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

## ---------- Load tAge Predictions ----------
tage_data <- read.csv("data/tage/Predicted_tAges_Gaucher_dataset.csv")
colnames(tage_data) <- gsub("^\uFEFF", "", colnames(tage_data))

tage_data$Group <- factor(
  ifelse(tage_data$Group == "Control", "Control",
  ifelse(tage_data$Group == "Gaucher", "GD", "GD+MA-5")),
  levels = c("Control", "GD", "GD+MA-5"))

group_colors <- c("Control" = "#6BAED6", "GD" = "#FC9272", "GD+MA-5" = "#74C476")

celltype_order <- c("microglia", "neuron", "astrocyte",
                    "oligodendrocyte", "opc", "endothelial cell")
tage_data$Cell.type <- factor(tolower(tage_data$Cell.type), levels = celltype_order)

## ---------- Dotplot Function ----------
create_tage_dotplot <- function(data, value_col, title) {
  df <- data.frame(
    Celltype = data$Cell.type, Group = data$Group, Value = data[[value_col]]
  ) %>% filter(!is.na(Value))

  plots <- list()
  for (ct in levels(df$Celltype)) {
    ct_data <- df %>% filter(Celltype == ct)
    if (nrow(ct_data) == 0) next
    p <- ggplot(ct_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6) +
      geom_jitter(aes(color = Group), width = 0.15, size = 2.5, alpha = 0.8) +
      scale_fill_manual(values = group_colors) +
      scale_color_manual(values = group_colors) +
      labs(title = tools::toTitleCase(ct), y = title, x = NULL) +
      theme_bw(base_family = "Arial") +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    plots[[ct]] <- p
  }
  wrap_plots(plots, nrow = 1)
}

create_tage_dotplot(tage_data, "tAge_chronological", "Chronological tAge")
create_tage_dotplot(tage_data, "tAge_mortality",     "Mortality tAge")

## ---------- Gene Contribution Barplots (Figure 7d) ----------
gene_contrib_gd  <- read.csv("data/tage/Gene_contribution_mortality_microglia_Gaucher.csv",
                              row.names = 1)
gene_contrib_ma5 <- read.csv("data/tage/Gene_contribution_mortality_microglia_GaucherMA5.csv",
                              row.names = 1)

prep_barplot_data <- function(df, n_top = 24) {
  df$abs_slope <- abs(df$Slope.Clock)
  df <- df[order(-df$abs_slope), ][1:n_top, ]
  df[order(df$Slope.Clock), ]
}

color_up   <- "#E74C3C"
color_down <- "#5DADE2"

create_contribution_barplot <- function(df, title) {
  df$gene_f <- factor(df$Gene.symbol, levels = df$Gene.symbol)
  ggplot(df, aes(x = gene_f, y = Slope.Clock, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = Slope.Clock - SE.Clock,
                      ymax = Slope.Clock + SE.Clock), width = 0.2) +
    scale_fill_manual(values = c("Upregulated" = color_up, "Downregulated" = color_down)) +
    coord_flip() +
    labs(title = title, x = NULL, y = expression("Effect on tAge, " * log[10] * "HR")) +
    theme_bw(base_family = "Arial") +
    theme(axis.text.y = element_text(size = 10))
}

create_contribution_barplot(prep_barplot_data(gene_contrib_gd),       "GD vs Control")
create_contribution_barplot(prep_barplot_data(gene_contrib_ma5, 25),  "GD+MA-5 vs GD")
