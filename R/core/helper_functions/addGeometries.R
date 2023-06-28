#'#' Add centroid and hexagon geometries to SpatialFeatureExperiment object
#'
#' @name addGeometries
#'
#' @description
#' A function to add centroid and hexagon geometries in the \code{colGeometries}
#' slot of a SpatialFeatureExperiment (SFE) object.
#'
#' @param sfe The SpatialFeatureExperiment object.
#'
#' @param samples A character vector specifying one or more directories, each
#' corresponding to a 10x Genomics Visium sample (see Details). If provided,
#' the names will be used as sample identifiers.
#'
#' @param sample_id A character string specifying unique sample identifiers,
#' one for each directory specified via samples. This parameter is ignored if
#' the names of samples are not null (!is.null(names(samples))).
#'
#' @param res The desired resolution. Can take one of "lowres", "hires",
#' or "fullres".
#'
#' @details
#' This function adds centroid and hexagon geometries to the
#' \code{colGeometries} slot of the SpatialFeatureExperiment object. It
#' calculates the spot diameter and adds spot centroids and hexagon geometries
#' accordingly.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso \code{\link{add.spotCntd}}, \code{\link{spot.diameter}},
#' \code{\link{add.spotHex}}, [read10xVisiumSFE()]
#'
#' @examples
#' # Load the SpatialFeatureExperiment object
#' data(sfe)
#'
#' # Set the file path to the sample folder.
#' # These are the folder paths you used when generating the SFE object using
#' # `read10xVisiumSFE` function.
#' sampleDirs <- "path/to/folder/"
#'
#' # Add geometries to the object
#' \dontrun{
#' sfe_object <- addGeometries(sfe = sfe,
#'                             samples = sampleDirs,
#'                             sample_id = c("sampleA", "sampleB"),
#'                             res = "hires")
#' }
#'
#' @export
addGeometries <- function(sfe,
                          samples,
                          sample_id,
                          res = c("lowres", "hires", "fullres")) {

  res <- match.arg(res)

  ## Add Centroids
  sfe <- add.spotCntd(sfe)

  ## Get/ calculate spot diameter
  sfe <- spot.diameter(sfe = sfe,
                       samples = samples,
                       sample_id = sample_id,
                       res = res)

  ## Add Hexagons
  sfe <- add.spotHex(sfe = sfe,
                     sample_id = sample_id,
                     res = res)

  return(sfe)
}
