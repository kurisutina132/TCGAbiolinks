\name{plot_pca}
\alias{plot_pca}
\title{PCA plot from a DESeq2 VST object}
\usage{
plot_pca(vsd, group_col = "sample_type")
}
\arguments{
\item{vsd}{A \code{DESeqTransform} object.}
\item{group_col}{Metadata column used for coloring samples.}
}
\value{A \code{ggplot} object.}
\description{Creates a PCA scatter plot from VST-transformed expression data.}
