% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_plot.R
\name{plotQC_spotsAnnotation}
\alias{plotQC_spotsAnnotation}
\title{Plot Quality Control of Spots with Annotation in a SpatialFeatureExperiment}
\usage{
plotQC_spotsAnnotation(
  m_sfe,
  sample_id = NULL,
  type = c("spot", "hex", "cntd"),
  fill_args = list(),
  colours = NULL,
  ...
)
}
\arguments{
\item{m_sfe}{A SpatialFeatureExperiment or MetaSFE object.}

\item{sample_id}{Character string, TRUE, or NULL specifying sample/image
identifier(s). TRUE is equivalent to all samples/images, and NULL specifies
the first available entry.}

\item{type}{Character vector specifying the spot type, either "spot", "cntd
or "hex".}

\item{fill_args}{List of arguments to customize the fill scale. Default is
an empty list.}

\item{colours}{Vector of colours for annotation. Default is NULL, which uses
default colours.}

\item{...}{Additional arguments to be passed to the underlying plotting
functions.}
}
\value{
A ggplot object representing the quality control plot with
annotation.
}
\description{
This function plots the quality control of spots in a
SpatialFeatureExperiment, including annotation information.
}
\details{
This function plots the quality control of spots in a
SpatialFeatureExperiment, including annotation information.

Users can specify the spot type (spots, centroids, or hexagons) and
customise the fill scale for annotation.

The colours argument allows users to customize annotation colours, and
additional arguments (...) can be passed to the underlying plotting
functions.
}
\examples{
\dontrun{
# Example usage:
data(sfe)
plotQC_spotsAnnotation(sfe, type = "spot", fill_args = list(),
                       sample_id = TRUE, colours = NULL)
}

}
\seealso{
Explore the related internal functions: \code{\link{.int_getSmplIDs}},
\code{\link{.int_dataToPlotAnnot}}, and \code{\link{getColours}}.
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{annotation}
\keyword{ggplot2}
\keyword{plotting}
\keyword{quality-control}
\keyword{spatial}
