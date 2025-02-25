% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalisation.R
\name{.int_msfeCompSizeFact}
\alias{.int_msfeCompSizeFact}
\title{INTERNAL: Compute Library Size Factors in a SpatialFeatureExperiment}
\usage{
.int_msfeCompSizeFact(sfe, ...)
}
\arguments{
\item{sfe}{A SpatialFeatureExperiment object.}

\item{...}{other arguments to be passed to \code{scater::computeLibraryFactors}}
}
\value{
A logical vector indicating the subset of features based on the
criteria.
}
\description{
INTERNAL: Compute Library Size Factors in a SpatialFeatureExperiment
}
