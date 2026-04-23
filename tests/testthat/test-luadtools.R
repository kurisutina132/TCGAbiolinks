test_that(".clean_ensembl_ids strips version suffix", {
  expect_equal(
    luadtools:::.clean_ensembl_ids(c("ENSG000001.12", "ENSG000002")),
    c("ENSG000001", "ENSG000002")
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
