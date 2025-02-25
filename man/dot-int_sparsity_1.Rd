% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_utils.R
\name{dot-int_sparsity_1}
\alias{dot-int_sparsity_1}
\alias{.int_sparsity_1}
\title{Internal: calculate feature sparsity}
\usage{
.int_sparsity_1(sfe, assay, sampleNo, .sample_id)
}
\arguments{
\item{sfe}{a SpatialFeaturesExperiment objects.}

\item{assay}{the name of the assay to use.}

\item{sampleNo}{used to call the sparsity also over aall datasets when
multiple samples exist.Used only when MARGIN = 1 at \code{addPerGeneQC}
function.}

\item{.sample_id}{the sample id to run the calculations for. Used only when
MARGIN = 1 at \code{addPerGeneQC} function.}
}
\description{
IMPORTANT: this function is not exported to be used alone.
A function to calculate the sparsity of the dataset and get the
sparsity of each gene.
}
\keyword{internal}
