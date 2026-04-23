\name{plot_topvar_heatmap}
\alias{plot_topvar_heatmap}
\title{Heatmap of top variable genes from a VST matrix}
\usage{
plot_topvar_heatmap(vsd, group, top_n = 20, ...)
}
\arguments{
\item{vsd}{A \code{DESeqTransform} object.}
\item{group}{Group annotation vector of length \code{ncol(vsd)}.}
\item{top_n}{Number of top variable genes to display.}
\item{...}{Additional arguments forwarded to \code{pheatmap::pheatmap()}.}
}
\value{The \code{pheatmap} return object.}
\description{Selects top variable genes, scales by row, and plots a sample-annotated heatmap.}
