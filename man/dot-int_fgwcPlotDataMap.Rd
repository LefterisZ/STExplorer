% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgwc_plot.R
\name{.int_fgwcPlotDataMap}
\alias{.int_fgwcPlotDataMap}
\title{Internal Function: Generate Data for FGWC Plot}
\usage{
.int_fgwcPlotDataMap(fgwc, sfe, mode = c("single", "multi"))
}
\arguments{
\item{fgwc}{An object containing FGWC results, typically obtained from the
\code{fgwcSTE} function.}

\item{sfe}{A SpatialFeatureExperiment object containing spatial coordinates
and features.}

\item{mode}{Character, indicating the mode for generating plot data. Options
are "single" (default) for a single-cluster mode or "multi" for a
multi-cluster mode.}
}
\value{
A data frame containing information for generating FGWC plots. In
single-cluster mode, the data frame includes cluster information for each
spot. In multi-cluster mode, the entire data frame is returned.
}
\description{
This function generates data suitable for plotting FFGWC results.
}
\details{
The function generates data suitable for plotting FGWC results. In
single-cluster mode, the data frame includes cluster information for each
spot. In multi-cluster mode, the entire data frame is returned. The mode is
specified using the 'mode' argument.
}
\seealso{
\code{\link[SpatialFeatureExperiment]{colGeometry}},
\code{\link[SpatialFeatureExperiment]{colData}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{FGWC,}
\keyword{clustering,}
\keyword{coordinates}
\keyword{spatial}
