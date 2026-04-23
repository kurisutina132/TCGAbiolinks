# luadtools

`luadtools` turns the TCGA-LUAD workflow into a reusable R package for interview-ready demonstration: data download/preparation, DESeq2 differential expression, and standard LUAD visual outputs.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kurisutina132/TCGAbiolinks")
```

## Minimal example

```r
library(luadtools)

# Network-heavy data retrieval
se <- tcga_luad_download_prepare(
  sample_types = c("Primary Tumor", "Solid Tissue Normal"),
  project = "TCGA-LUAD"
)

out <- deseq2_tumor_vs_normal(se)

p_pca <- plot_pca(out$vsd, group_col = "sample_type")
p_volcano <- plot_volcano(out$res)
p_heatmap <- plot_topvar_heatmap(out$vsd, group = out$coldata$sample_type, top_n = 20)

kegg <- kegg_enrichment_from_deseq(out$res, direction = "up")
```

## Outputs

- `tcga_luad_download_prepare()`: `SummarizedExperiment` LUAD object
- `deseq2_tumor_vs_normal()`: list with `dds`, `res`, `vsd`, `coldata`
- `plot_pca()`: PCA `ggplot`
- `plot_volcano()`: volcano `ggplot`
- `plot_topvar_heatmap()`: `pheatmap` object
- `kegg_enrichment_from_deseq()`: list with `enrichment`, `dotplot`, `barplot`, `genes`

## Troubleshooting

- **ENSEMBL version suffix**: use IDs without suffix (e.g., `ENSG...` not `ENSG....12`) for mapping; package helper removes suffix when needed.
- **Factor levels**: set normal as reference (`Solid Tissue Normal`) before DE analysis to interpret tumor-vs-normal correctly.
- **`padj` NA values**: common in low-information genes; filter using `!is.na(padj)` before significance summaries.
- **Download issues**: GDC servers can be rate-limited or temporarily unavailable; retry later and keep stable internet.
- **Reproducibility**: use `renv::init()` and lock package versions before sharing results.

## Reference script

Original workflow is preserved at:

`inst/scripts/TCGAbiolinks_LUAD_workflow.Rmd`
