% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_plot.R
\name{.int_getLocalResName}
\alias{.int_getLocalResName}
\title{Internal Function: Get localResults Name for Spatial Autocorrelation}
\usage{
.int_getLocalResName(statistic, test)
}
\arguments{
\item{statistic}{A character string specifying the type of spatial
autocorrelation used ("moran", "geary", "getis").}

\item{test}{A character string specifying the statistical test method used
("z-score" or "permutation").}
}
\value{
A character string representing the name of the local results.
}
\description{
Internal Function: Get localResults Name for Spatial Autocorrelation
}
\author{
Eleftherios (Lefteris) Zormpas
}
