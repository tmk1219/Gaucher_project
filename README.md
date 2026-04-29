# Gaucher_project

Analysis code for:

**Mitochondria-targeted restoration of mtDNA integrity suppresses cGAS-STING inflammation and rejuvenates Gaucher disease phenotypes**

## Data Availability

RNA-seq data are deposited in NCBI GEO under SuperSeries [GSE325369](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325369).

## Scripts

| Script | Description | Figure |
|--------|-------------|--------|
| `01_snRNAseq_preprocessing.R` | QC, integration, clustering, cell type annotation | Fig 3b-d |
| `02_pathway_activity_scoring.R` | MSigDB Hallmark pathway scoring per cell | Fig 3e |
| `03_signature_scores.R` | DAM, IFN, Activation, NDMG, SASP scores | Fig 4 |
| `04_microglia_subclustering.R` | Microglia re-clustering and proportion analysis | Fig 5a-c |
| `05_microglia_marker_enrichment.R` | Cluster markers and Enrichr pathway analysis | Fig 5d-g |
| `06_scafe_cre_quantification.R` | SCAFE downstream tCRE activity quantification | Fig 6a |
| `07_scafe_genome_tracks.R` | Gviz genome track visualization | Fig 6c-i |
| `08_bulk_rnaseq_downstream.R` | PCA, rescue gene identification, Enrichr | Fig 8p,r, 9d-f |
| `09_bulk_gsea.R` | fgsea pathway enrichment on bulk RNA-seq | Fig 8q, 9g |
| `10_bulk_volcano_plots.R` | Volcano plot generation | Supplementary |
| `11_tage_analysis.R` | Transcriptomic age metacell and clock pipeline | Fig 7a,c |
| `12_tage_visualization.R` | tAge dotplots and gene contribution barplots | Fig 7a,c,d |
| `13_correlation_heatmap.R` | Spearman correlation heatmap | Fig 7b |
| `14_gsea_heatmap.R` | GSEA NES heatmap visualization | Fig 7b, 8q, 9g |

## Dependencies

- R (v4.4.2): Seurat (v5.3.1), MAST (v1.22.0), fgsea (v1.34.2), enrichR (v3.2), Gviz, rtracklayer, metafor
- Bulk RNA-seq differential expression was performed using CLC Genomics Workbench
- Transcriptomic age (tAge) analysis uses the [tAge framework](https://github.com/Gladyshev-Lab/tAge) (Tyshkovskiy et al., 2024)
- SCAFE processing was performed using Docker
