\name{deseq2_tumor_vs_normal}
\alias{deseq2_tumor_vs_normal}
\title{Differential expression: tumor vs normal with DESeq2}
\usage{
deseq2_tumor_vs_normal(
  se,
  sample_col = "sample_type",
  tumor_label = "Primary Tumor",
  normal_label = "Solid Tissue Normal",
  min_total_counts = 10,
  shrink_type = "apeglm"
)
}
\arguments{
\item{se}{A \code{SummarizedExperiment} with counts and sample metadata.}
\item{sample_col}{Metadata column indicating sample group.}
\item{tumor_label}{Tumor level label in \code{sample_col}.}
\item{normal_label}{Normal/reference level label in \code{sample_col}.}
\item{min_total_counts}{Minimum total count per gene to keep.}
\item{shrink_type}{Shrinkage method passed to \code{DESeq2::lfcShrink()}.}
}
\value{A named list with \code{dds}, \code{res}, \code{vsd}, and \code{coldata}.}
\description{Builds a DESeq2 dataset, sets normal as reference level, runs DESeq once, shrinks LFC, and computes VST output.}
\examples{
\dontrun{
se <- tcga_luad_download_prepare()
out <- deseq2_tumor_vs_normal(se)
}
}
