\name{kegg_enrichment_from_deseq}
\alias{kegg_enrichment_from_deseq}
\title{KEGG enrichment from DESeq2 results}
\usage{
kegg_enrichment_from_deseq(
  res,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  direction = c("all", "up", "down")
)
}
\arguments{
\item{res}{Result object/data frame containing \code{log2FoldChange} and \code{padj}.}
\item{padj_cutoff}{Adjusted p-value cutoff.}
\item{lfc_cutoff}{Absolute log2 fold-change cutoff.}
\item{direction}{One of \code{"all"}, \code{"up"}, or \code{"down"}.}
}
\value{A list with \code{enrichment}, \code{dotplot}, \code{barplot}, and \code{genes}.}
\description{Filters DE genes, removes Ensembl version suffixes, maps to ENTREZ IDs, and runs KEGG enrichment with dotplot/barplot outputs.}
