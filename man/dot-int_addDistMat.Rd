% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.int_addDistMat}
\alias{.int_addDistMat}
\title{Internal Function: Add Distance Matrix to SpatialFeatureExperiment}
\usage{
.int_addDistMat(sfe, p, ...)
}
\arguments{
\item{sfe}{A SpatialFeatureExperiment object containing spatial coordinates.}

\item{p}{Numeric, the power parameter for the distance metric (default is 2
for Euclidean distance).}

\item{...}{Additional arguments to be passed to the gw.dist function.}
}
\value{
A modified SpatialFeatureExperiment object with the distance matrix
added to its metadata.
}
\description{
This function calculates the distance matrix for spatial coordinates and
adds it to the metadata of a SpatialFeatureExperiment object.
}
\details{
The function calculates the distance matrix for the spatial
coordinates in the SpatialFeatureExperiment object using the specified
distance metric (Euclidean, Manhattan, or Minkowski). The resulting distance
matrix is stored in the metadata of the SpatialFeatureExperiment object.
}
\seealso{
\code{\link[GWmodel]{gw.dist}},
\code{\link[SpatialFeatureExperiment]{spatialCoords}},
\code{\link[SpatialFeatureExperiment]{colData}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{SpatialFeatureExperiment}
\keyword{coordinates,}
\keyword{distance}
\keyword{matrix,}
\keyword{spatial}
