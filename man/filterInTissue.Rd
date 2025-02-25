% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_filterInTissue.R
\name{filterInTissue}
\alias{filterInTissue}
\title{Filter In-Tissue Spots in Spatial Transcriptomics Data}
\usage{
filterInTissue(sfe, sample_id = TRUE)
}
\arguments{
\item{sfe}{A \code{SpatialFeatureExperiment} object.}

\item{sample_id}{A character string, TRUE, or NULL specifying sample/image
identifier(s).}
}
\value{
A \code{SpatialFeatureExperiment} object with in-tissue spots
filtered.
}
\description{
This function filters in-tissue spots in spatial transcriptomics data,
providing users with the flexibility to specify sample/image identifiers.
The sample_id parameter accepts TRUE (equivalent to all samples/images) or
NULL (specifying the first available entry).
}
\details{
The function internally utilises the .int_filterInTissue function to filter
in-tissue spots for each specified sample. The results are then merged into a
single SpatialFeatureExperiment object using the Reduce function.

Additionally, the function ensures data integrity by cleaning up duplicate
entries in the metadata slot with the .int_cleanMetaData function.

This filtering step is crucial for focusing analyses on relevant in-tissue
spots in spatial transcriptomics experiments.
}
\examples{
\dontrun{
sfe <- filterInTissue(sfe)
}


}
\seealso{
\code{\link{.int_filterInTissue}}, \code{\link{.int_cleanMetaData}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{data}
\keyword{filtering}
\keyword{in-tissue}
\keyword{spatial}
\keyword{spots,}
\keyword{transcriptomics,}
