test_that(".clean_ensembl_ids strips version suffix", {
  expect_equal(
    luadtools:::.clean_ensembl_ids(c("ENSG000001.12", "ENSG000002")),
    c("ENSG000001", "ENSG000002")
  )
})

test_that("deseq2_tumor_vs_normal validates input type", {
  expect_error(
    deseq2_tumor_vs_normal(iris),
    "`se` must be a SummarizedExperiment."
  )
})

test_that("plot_volcano returns a ggplot object", {
  mock_res <- data.frame(
    log2FoldChange = c(2, -1.5, 0.1),
    padj = c(0.001, 0.01, 0.9),
    row.names = c("ENSG000001.12", "ENSG000002.3", "ENSG000003.7")
  )

  p <- plot_volcano(mock_res, title = "Test", subtitle = "Mock")
  expect_s3_class(p, "ggplot")
})

test_that("plot_topvar_heatmap returns pheatmap object with matrix-like assay", {
  mat <- matrix(
    c(3, 2, 4, 9, 8, 10, 5, 6, 7, 1, 2, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("ENSG000001.12", "ENSG000002.3", "ENSG000003.7"),
      c("S1", "S2", "S3", "S4")
    )
  )
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  ph <- plot_topvar_heatmap(se, group = c("Tumor", "Tumor", "Normal", "Normal"), top_n = 2)
  expect_s3_class(ph, "pheatmap")
})

test_that("kegg_enrichment_from_deseq validates required columns", {
  expect_error(
    kegg_enrichment_from_deseq(data.frame(x = 1)),
    "`res` must contain `log2FoldChange` and `padj` columns."
  )
})
