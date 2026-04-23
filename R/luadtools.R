#' Download and prepare TCGA-LUAD RNA-seq counts
#'
#' Wraps `TCGAbiolinks::GDCquery()`, `TCGAbiolinks::GDCdownload()`, and
#' `TCGAbiolinks::GDCprepare()` for LUAD RNA-seq counts.
#'
#' @param sample_types Character vector of TCGA sample types.
#' @param project TCGA project identifier.
#'
#' @return A `SummarizedExperiment` object.
#' @export
#'
#' @examples
#' \dontrun{
#' se <- tcga_luad_download_prepare()
#' }
tcga_luad_download_prepare <- function(
  sample_types = c("Primary Tumor", "Solid Tissue Normal"),
  project = "TCGA-LUAD"
) {
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = sample_types
  )

  TCGAbiolinks::GDCdownload(query)
  TCGAbiolinks::GDCprepare(query)
}

#' Differential expression: tumor vs normal with DESeq2
#'
#' Builds a DESeq2 dataset from a `SummarizedExperiment`, enforces a
#' normal-reference factor level, runs DESeq once, applies LFC shrinkage,
#' and computes a variance-stabilized transform.
#'
#' @param se A `SummarizedExperiment` with counts and sample metadata.
#' @param sample_col Metadata column indicating sample group.
#' @param tumor_label Tumor level label in `sample_col`.
#' @param normal_label Normal/reference level label in `sample_col`.
#' @param min_total_counts Minimum total count per gene to keep.
#' @param shrink_type Shrinkage method passed to `DESeq2::lfcShrink()`.
#'
#' @return A named list with `dds`, `res`, `vsd`, and `coldata`.
#' @export
#'
#' @examples
#' \dontrun{
#' se <- tcga_luad_download_prepare()
#' out <- deseq2_tumor_vs_normal(se)
#' }
deseq2_tumor_vs_normal <- function(
  se,
  sample_col = "sample_type",
  tumor_label = "Primary Tumor",
  normal_label = "Solid Tissue Normal",
  min_total_counts = 10,
  shrink_type = "apeglm"
) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }

  counts <- SummarizedExperiment::assay(se)
  coldata <- as.data.frame(SummarizedExperiment::colData(se))

  if (!sample_col %in% names(coldata)) {
    stop("`sample_col` not found in colData(se).", call. = FALSE)
  }

  groups <- as.character(coldata[[sample_col]])
  if (!all(c(tumor_label, normal_label) %in% unique(groups))) {
    stop("`tumor_label` and/or `normal_label` not found in sample metadata.", call. = FALSE)
  }

  coldata[[sample_col]] <- factor(groups, levels = c(normal_label, tumor_label))

  design_formula <- stats::as.formula(paste0("~", sample_col))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = coldata,
    design = design_formula
  )

  keep <- base::rowSums(DESeq2::counts(dds)) >= min_total_counts
  dds <- dds[keep, ]

  dds <- DESeq2::DESeq(dds)

  result_names <- DESeq2::resultsNames(dds)
  tumor_tag <- make.names(tumor_label)
  normal_tag <- make.names(normal_label)
  contrast_coef <- result_names[
    grepl("_vs_", result_names) &
      grepl(tumor_tag, result_names, fixed = TRUE) &
      grepl(normal_tag, result_names, fixed = TRUE)
  ]

  if (length(contrast_coef) == 0L) {
    stop("Could not identify tumor-vs-normal coefficient in DESeq2 resultsNames.", call. = FALSE)
  }

  res <- DESeq2::lfcShrink(dds, coef = contrast_coef[[1]], type = shrink_type)
  vsd <- DESeq2::vst(dds, blind = FALSE)

  list(
    dds = dds,
    res = res,
    vsd = vsd,
    coldata = as.data.frame(SummarizedExperiment::colData(dds))
  )
}

#' PCA plot from a DESeq2 VST object
#'
#' @param vsd A `DESeqTransform` object.
#' @param group_col Metadata column used for coloring samples.
#'
#' @return A `ggplot` object.
#' @export
plot_pca <- function(vsd, group_col = "sample_type") {
  pca_data <- DESeq2::plotPCA(vsd, intgroup = group_col, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  pca_data$group <- pca_data[[group_col]]

  ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, color = group)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggplot2::labs(color = group_col) +
    ggplot2::theme_bw(base_size = 12)
}

#' Volcano plot from DESeq2 results
#'
#' @param res Result object/data frame containing `log2FoldChange` and `padj`.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#'
#' @return A `ggplot` object produced by EnhancedVolcano.
#' @export
plot_volcano <- function(
  res,
  title = "TCGA-LUAD: Tumor vs Normal",
  subtitle = "Transcriptomic biomarkers"
) {
  res_df <- as.data.frame(res)

  EnhancedVolcano::EnhancedVolcano(
    res_df,
    lab = rownames(res_df),
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    title = title,
    subtitle = subtitle,
    labSize = 3,
    legendPosition = "right"
  )
}

#' Heatmap of top variable genes from a VST matrix
#'
#' @param vsd A `DESeqTransform` object.
#' @param group Group annotation vector of length `ncol(vsd)`.
#' @param top_n Number of top variable genes to display.
#' @param ... Additional arguments forwarded to `pheatmap::pheatmap()`.
#'
#' @return The `pheatmap` return object.
#' @export
plot_topvar_heatmap <- function(vsd, group, top_n = 20, ...) {
  mat <- SummarizedExperiment::assay(vsd)
  if (length(group) != ncol(mat)) {
    stop("`group` must have length equal to ncol(vsd).", call. = FALSE)
  }
  vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  top_idx <- utils::head(order(vars, decreasing = TRUE), top_n)

  mat_top <- mat[top_idx, , drop = FALSE]
  mat_top_scaled <- t(scale(t(mat_top)))
  mat_top_scaled[is.na(mat_top_scaled)] <- 0

  annotation <- data.frame(group = group)
  rownames(annotation) <- colnames(mat_top_scaled)

  pheatmap::pheatmap(
    mat_top_scaled,
    annotation_col = annotation,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    ...
  )
}

#' KEGG enrichment from DESeq2 results
#'
#' @param res Result object/data frame containing `log2FoldChange` and `padj`.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param lfc_cutoff Absolute log2 fold-change cutoff.
#' @param direction One of `"all"`, `"up"`, or `"down"`.
#'
#' @return A list with `enrichment`, `dotplot`, `barplot`, and `genes`.
#' @export
kegg_enrichment_from_deseq <- function(
  res,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  direction = c("all", "up", "down")
) {
  direction <- match.arg(direction)

  res_df <- as.data.frame(res)
  if (!all(c("log2FoldChange", "padj") %in% names(res_df))) {
    stop("`res` must contain `log2FoldChange` and `padj` columns.", call. = FALSE)
  }

  res_df$gene <- rownames(res_df)

  keep <- !is.na(res_df$padj) &
    (res_df$padj <= padj_cutoff) &
    (abs(res_df$log2FoldChange) >= lfc_cutoff)

  if (direction == "up") {
    keep <- keep & res_df$log2FoldChange > 0
  } else if (direction == "down") {
    keep <- keep & res_df$log2FoldChange < 0
  }

  genes <- .clean_ensembl_ids(res_df$gene[keep])
  genes <- unique(stats::na.omit(genes))

  if (length(genes) == 0L) {
    return(list(enrichment = NULL, dotplot = NULL, barplot = NULL, genes = character(0)))
  }

  entrez <- clusterProfiler::bitr(
    genes,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  if (is.null(entrez) || nrow(entrez) == 0L) {
    return(list(enrichment = NULL, dotplot = NULL, barplot = NULL, genes = genes))
  }

  enrichment <- clusterProfiler::enrichKEGG(
    gene = entrez$ENTREZID,
    organism = "hsa",
    pvalueCutoff = padj_cutoff
  )

  if (is.null(enrichment) || nrow(as.data.frame(enrichment)) == 0L) {
    return(list(enrichment = enrichment, dotplot = NULL, barplot = NULL, genes = genes))
  }

  dot <- clusterProfiler::dotplot(enrichment, showCategory = 15) +
    ggplot2::ggtitle(paste("KEGG enrichment -", direction, "genes"))

  bar <- clusterProfiler::barplot(
    enrichment,
    showCategory = 15,
    title = paste("KEGG enrichment -", direction, "genes")
  )

  list(enrichment = enrichment, dotplot = dot, barplot = bar, genes = genes)
}

.clean_ensembl_ids <- function(ids) {
  gsub("\\..*", "", ids)
}
