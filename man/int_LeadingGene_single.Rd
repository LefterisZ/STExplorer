% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwpca_LeadingGene.R
\name{int_LeadingGene_single}
\alias{int_LeadingGene_single}
\title{Find the single leading gene}
\usage{
int_LeadingGene_single(gwpca, pc_no, sfe, sort)
}
\arguments{
\item{gwpca}{a list of class \code{gwpca}.}

\item{pc_no}{a numeric value of the principal component number for which
you want to find the leading genes.}

\item{sfe}{a \code{SpatialFeatureExperiment} object.}

\item{sort}{A character indicating whether the absolute or the original
leading score will be used for sorting. Possible values are "abs" (for
absolute) or "original"}
}
\description{
A function to find the SINGLE leading gene per principal component in a
location.
}
\keyword{internal}
