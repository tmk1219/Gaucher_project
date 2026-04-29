## =========================================================
##  SCAFE Downstream: CRE Activity Quantification
##  Figure 6a
##
##  SCAFE was run via Docker for TSS cluster identification.
##  This script quantifies tCRE activity from SCAFE output.
##
##  Input:  bigWig files (TPM-normalized SCAFE output)
##  Output: tCRE activity matrices
## =========================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

## ---------- Load SCAFE Output ----------
bw_ctrl <- import("data/scafe/Control_microglia.bw")
bw_gd   <- import("data/scafe/GD_microglia.bw")
bw_ma5  <- import("data/scafe/MA5_microglia.bw")

## ---------- Load tCRE Classification Parameters ----------
tcre_params <- read.delim("config/tcre_classification_parameters.txt", header = TRUE)

## ---------- Load FANTOM5 Annotations ----------
fantom5_enhancers <- import.bed(tcre_params$value[tcre_params$parameter == "fantom5_enhancer_bed"])
fantom5_promoters <- import.bed(tcre_params$value[tcre_params$parameter == "fantom5_promoter_bed"])

## ---------- Quantify CRE Activity ----------
quantify_cre <- function(bw, regions, mode = "promoter") {
  scores <- rep(0, length(regions))
  for (i in seq_along(regions)) {
    region <- regions[i]
    if (mode == "promoter") {
      # Max signal for sharp TSS peaks
      ov <- subsetByOverlaps(bw, region)
      scores[i] <- ifelse(length(ov) > 0, max(ov$score), 0)
    } else {
      # Mean signal within +/-1kb for enhancers
      expanded <- resize(region, width = 2000, fix = "center")
      ov <- subsetByOverlaps(bw, expanded)
      scores[i] <- ifelse(length(ov) > 0, mean(ov$score), 0)
    }
  }
  scores
}

min_tpm <- as.numeric(tcre_params$value[tcre_params$parameter == "min_tpm_threshold"])

# Promoter activity
prom_ctrl <- quantify_cre(bw_ctrl, fantom5_promoters, mode = "promoter")
prom_gd   <- quantify_cre(bw_gd,   fantom5_promoters, mode = "promoter")
prom_ma5  <- quantify_cre(bw_ma5,  fantom5_promoters, mode = "promoter")

# Enhancer activity
enh_ctrl <- quantify_cre(bw_ctrl, fantom5_enhancers, mode = "enhancer")
enh_gd   <- quantify_cre(bw_gd,   fantom5_enhancers, mode = "enhancer")
enh_ma5  <- quantify_cre(bw_ma5,  fantom5_enhancers, mode = "enhancer")

# Filter by minimum TPM
prom_mean <- (prom_ctrl + prom_gd + prom_ma5) / 3
enh_mean  <- (enh_ctrl + enh_gd + enh_ma5) / 3

message(sprintf("Promoter tCREs retained: %d", sum(prom_mean >= min_tpm)))
message(sprintf("Enhancer tCREs retained: %d", sum(enh_mean >= min_tpm)))
