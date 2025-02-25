% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.int_getSmplIDs}
\alias{.int_getSmplIDs}
\title{Internal Function: .int_getSmplIDs}
\usage{
.int_getSmplIDs(sfe, sample_id = NULL)
}
\arguments{
\item{sfe}{A SpatialFeatureExperiment object.}

\item{sample_id}{A character string, \code{TRUE}, or \code{NULL} specifying
sample/image identifier(s). If \code{TRUE}, all samples/images are
considered. If \code{NULL}, the first available entry is considered.}
}
\value{
Returns a character vector of sample IDs for plotting.
}
\description{
Description: Fetches sample IDs.
}
\details{
This function fetches sample IDs based on the input parameters,
providing flexibility for customizing the sample IDs for plotting.
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{internal}
