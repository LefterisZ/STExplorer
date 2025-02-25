% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_plot.R
\name{.int_plotSpot}
\alias{.int_plotSpot}
\title{Internal Function: .int_plotSpot}
\usage{
.int_plotSpot(sfe, in_tissue, n_samples, colours, ids)
}
\arguments{
\item{sfe}{A SpatialFeatureExperiment object.}

\item{in_tissue}{Logical. If TRUE, plot only in-tissue spots; if FALSE,
plot all spots.}

\item{n_samples}{Number of samples in the dataset.}

\item{colours}{Vector of colours for filling the spots. If NULL, default
colours are used.}

\item{ids}{Vector of sample IDs to include in the plot.}
}
\value{
Returns a ggplot2 plot object.
}
\description{
Generates a ggplot2 plot of spatial coordinates for spots with
optional colour customisation.
}
\details{
This function combines sample metadata with spatial coordinates and
creates a plot using ggplot2.
The plot can be customised based on the input parameters.
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{internal}
