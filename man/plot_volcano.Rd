\name{plot_volcano}
\alias{plot_volcano}
\title{Volcano plot from DESeq2 results}
\usage{
plot_volcano(
  res,
  title = "TCGA-LUAD: Tumor vs Normal",
  subtitle = "Transcriptomic biomarkers"
)
}
\arguments{
\item{res}{Result object/data frame containing \code{log2FoldChange} and \code{padj}.}
\item{title}{Plot title.}
\item{subtitle}{Plot subtitle.}
}
\value{A \code{ggplot} object produced by EnhancedVolcano.}
\description{Builds a standard volcano plot for LUAD DESeq2 results.}
