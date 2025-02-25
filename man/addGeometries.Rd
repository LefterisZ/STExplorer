% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_addGeometries.R
\name{addGeometries}
\alias{addGeometries}
\title{#' Add centroid and hexagon geometries to SpatialFeatureExperiment object}
\usage{
addGeometries(
  m_sfe,
  samples,
  sample_id,
  res = c("lowres", "hires", "fullres"),
  geoms = c("centroids", "hexagons", "both"),
  flipped = FALSE,
  barcodes = "all"
)
}
\arguments{
\item{m_sfe}{An object of class SpatialFeatureExperiment or
MetaSpatialFeatureExperiment.}

\item{samples}{A character vector specifying one or more directories, each
corresponding to a 10x Genomics Visium sample (see Details). If provided,
the names will be used as sample identifiers.}

\item{sample_id}{A character string specifying unique sample identifiers,
one for each directory specified via samples. This parameter is ignored if
the names of samples are not null (!is.null(names(samples))).}

\item{res}{The desired resolution. Can take one of "lowres", "hires",
or "fullres".}

\item{geoms}{A character string specifying the geometry types to be
generated. It has to be one of \code{centroids}, \code{hexagons}, or \code{both}. If left
empty, defaults to \code{both}. We suggest you leave it empty unless you are
analysing Curio Seeker (Slide-seq) data.}

\item{flipped}{Default is FALSE. This argument is important for 10X Visium
data. See details below.}

\item{barcodes}{A character string. Can take values either "all", "input",
or a user specified character vector with barcodes to be selected for
geometries generation. When you have an already pre-processed dataset, use
"input". Default is "all".}
}
\description{
A function to add centroid and hexagon geometries in the \code{colGeometries}
slot of a SpatialFeatureExperiment (SFE) object.
}
\details{
This function adds centroid and hexagon geometries to the
\code{colGeometries} slot of the SpatialFeatureExperiment object. It
calculates the spot diameter and adds spot centroids and hexagon geometries
accordingly.
About the \code{flipped} argument. Leave it as is and change it only if you see an
error like the one discussed in the vignette at the preprocessing step. 10X
Visium data sometimes are flipped. This means that the array layout and the
image pixel array are not on a same coordinate space but, most of the times,
the image is Y-flipped or rotated. If setting the flipped to \code{TRUE} doesn't
solve the occurring error, the please submit an issue at the GitHub repo of
STExplorer: https://github.com/LefterisZ/STExplorer/issues. The developers
will try to solve it for you.
}
\examples{
\dontrun{
# Load the SpatialFeatureExperiment object
data(sfe)

# Set the file path to the sample folder.
# These are the folder paths you used when generating the SFE object using
# `read10xVisiumSFE` function.
sampleDirs <- "path/to/folder/"

# Add geometries to the object
sfe_object <- addGeometries(sfe = sfe,
                            samples = sampleDirs,
                            sample_id = c("sampleA", "sampleB"),
                            res = "hires")
}

}
\seealso{
\code{\link{add.spotCntd}}, \code{\link{spot.diameter}},
\code{\link{add.spotHex}}, \code{\link[=read10xVisiumSFE]{read10xVisiumSFE()}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
