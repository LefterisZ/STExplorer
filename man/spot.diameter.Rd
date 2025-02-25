% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_utils.R
\name{spot.diameter}
\alias{spot.diameter}
\title{Calculate spot diameter}
\usage{
spot.diameter(sfe, samples, sample_id, res = c("lowres", "hires", "fullres"))
}
\arguments{
\item{sfe}{The SpatialFeaturesExperiment object.}

\item{samples}{A character vector specifying one or more directories, each
corresponding to a 10x Genomics Visium sample. If provided, the names will
be used as sample identifiers.}

\item{sample_id}{A character string specifying unique sample identifiers,
one for each directory specified via samples. This parameter is ignored if
names(samples) is not NULL.}

\item{res}{The resolution used to calculate the pixel XY coordinates for
each spot. It can take the values "lowres", "hires", or "fullres".}
}
\description{
IMPORTANT: this function is not exported to be used alone.
It is used internally to calculate the diameter, in pixels, of each Visium
slide spot using the scale factors. The scale factors can be found in the
"scalefactors_json.json" output from spaceranger. It is recommended to
provide the pathways to the folder where the "scalefactors_json.json" files
are located using the samples parameter generated using the file.path
function.
}
\details{
This function calculates the diameter, in pixels, of each Visium slide spot
using the scale factors provided in the "scalefactors_json.json" files.
The scale factors are imported from the specified directories (samples). The
res parameter specifies the resolution to use for the calculations:
"lowres", "hires", or "fullres". The calculated spot diameters are added to
the metadata of the sfe object.
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{internal}
