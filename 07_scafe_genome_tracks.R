## =========================================================
##  SCAFE Downstream: Genome Track Visualization
##  Figure 6c-i
##
##  Input:  bigWig files, genomic loci definitions
##  Output: Gviz genome track plots (PDF/PNG)
## =========================================================

suppressPackageStartupMessages({
  library(Gviz)
  library(rtracklayer)
  library(GenomicRanges)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

## ---------- Load bigWig Files ----------
bw_files <- list(
  Control   = "data/scafe/Control_microglia.bw",
  GD        = "data/scafe/GD_microglia.bw",
  `GD+MA-5` = "data/scafe/MA5_microglia.bw"
)

## ---------- Load Loci from Config ----------
loci <- read.delim("config/genomic_loci.txt", header = TRUE)

## ---------- Track Plotting Function ----------
plot_gene_tracks <- function(gene_symbol, chr, start, end, bw_files, txdb) {
  genome_track <- GenomeAxisTrack()
  gene_track <- GeneRegionTrack(txdb, chromosome = chr, start = start, end = end,
                                 name = gene_symbol, transcriptAnnotation = "symbol")

  data_tracks <- lapply(names(bw_files), function(cond) {
    DataTrack(range = bw_files[[cond]], chromosome = chr,
              from = start, to = end,
              name = cond, type = "h", col = "navy",
              genome = "mm10", ylim = c(0, NA))
  })

  plotTracks(c(genome_track, data_tracks, gene_track),
             from = start, to = end, chromosome = chr,
             main = gene_symbol, cex.main = 1.2)
}

## ---------- Generate All Tracks ----------
dir.create("output/genome_tracks", showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(loci))) {
  pdf(sprintf("output/genome_tracks/%s_Fig%s.pdf",
              loci$gene[i], loci$figure_panel[i]),
      width = 8, height = 6)
  plot_gene_tracks(loci$gene[i], loci$chr[i], loci$start[i], loci$end[i],
                   bw_files, txdb)
  dev.off()
}
