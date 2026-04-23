\name{tcga_luad_download_prepare}
\alias{tcga_luad_download_prepare}
\title{Download and prepare TCGA-LUAD RNA-seq counts}
\usage{
tcga_luad_download_prepare(
  sample_types = c("Primary Tumor", "Solid Tissue Normal"),
  project = "TCGA-LUAD"
)
}
\arguments{
\item{sample_types}{Character vector of TCGA sample types.}
\item{project}{TCGA project identifier.}
}
\value{A \code{SummarizedExperiment} object.}
\description{Wraps \code{GDCquery}, \code{GDCdownload}, and \code{GDCprepare} for LUAD RNA-seq counts.}
\examples{
\dontrun{tcga_luad_download_prepare()}
}
