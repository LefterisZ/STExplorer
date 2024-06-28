#'#' Add centroid and hexagon geometries to SpatialFeatureExperiment object
#'
#' @name addGeometries
#'
#' @description
#' A function to add centroid and hexagon geometries in the \code{colGeometries}
#' slot of a SpatialFeatureExperiment (SFE) object.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#'
#' @param samples A character vector specifying one or more directories, each
#' corresponding to a 10x Genomics Visium sample (see Details). If provided,
#' the names will be used as sample identifiers.
#'
#' @param sample_id A character string specifying unique sample identifiers,
#' one for each directory specified via samples. This parameter is ignored if
#' the names of samples are not null (!is.null(names(samples))).
#'
#' @param geoms A character string specifying the geometry types to be
#' generated. It has to be one of `centroids`, `hexagons`, or `both`. If left
#' empty, defaults to `both`. We suggest you leave it empty unless you are
#' analysing Curio Seeker (Slide-seq) data.
#'
#' @param res The desired resolution. Can take one of "lowres", "hires",
#' or "fullres".
#'
#' @param flipped Default is FALSE. This argument is important for 10X Visium
#' data. See details below.
#'
#' @param barcodes A character string. Can take values either "all", "input",
#' or a user specified character vector with barcodes to be selected for
#' geometries generation. When you have an already pre-processed dataset, use
#' "input". Default is "all".
#'
#' @details
#' This function adds centroid and hexagon geometries to the
#' \code{colGeometries} slot of the SpatialFeatureExperiment object. It
#' calculates the spot diameter and adds spot centroids and hexagon geometries
#' accordingly.
#' About the `flipped` argument. Leave it as is and change it only if you see an
#' error like the one discussed in the vignette at the preprocessing step. 10X
#' Visium data sometimes are flipped. This means that the array layout and the
#' image pixel array are not on a same coordinate space but, most of the times,
#' the image is Y-flipped or rotated. If setting the flipped to `TRUE` doesn't
#' solve the occurring error, the please submit an issue at the GitHub repo of
#' STExplorer: https://github.com/LefterisZ/STExplorer/issues. The developers
#' will try to solve it for you.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso \code{\link{add.spotCntd}}, \code{\link{spot.diameter}},
#' \code{\link{add.spotHex}}, [read10xVisiumSFE()]
#'
#' @examples
#' \dontrun{
#' # Load the SpatialFeatureExperiment object
#' data(sfe)
#'
#' # Set the file path to the sample folder.
#' # These are the folder paths you used when generating the SFE object using
#' # `read10xVisiumSFE` function.
#' sampleDirs <- "path/to/folder/"
#'
#' # Add geometries to the object
#' sfe_object <- addGeometries(sfe = sfe,
#'                             samples = sampleDirs,
#'                             sample_id = c("sampleA", "sampleB"),
#'                             res = "hires")
#' }
#'
#' @export
addGeometries <- function(m_sfe,
                          samples,
                          sample_id,
                          res = c("lowres", "hires", "fullres"),
                          geoms = c("centroids", "hexagons", "both"),
                          flipped = FALSE,
                          barcodes = "all") {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check arguments
  res <- match.arg(res)
  if (missing(geoms)) {
    geoms <- "both"
  }

  ## Determine which geoms to calculate based on the geoms argument
  if (geoms == "centroids") {
    doCentroids <- TRUE
    doHexagons <- FALSE
  } else if (geoms == "hexagons") {
    doCentroids <- FALSE
    doHexagons <- TRUE
  } else if (geoms == "both") {
    doCentroids <- TRUE
    doHexagons <- TRUE
  } else {
    stop("Invalid parameter value. Please use 'centroids',",
         " 'hexagons', or 'both'.")
  }

  ## Add Centroids
  if (doCentroids) {
    sfe <- add.spotCntd(sfe,
                        sample_id = sample_id)
  }

  if (doHexagons) {
    ## Get/ calculate spot diameter
    sfe <- spot.diameter(sfe = sfe,
                         samples = samples,
                         sample_id = sample_id,
                         res = res)

    ## Add Hexagons
    # NOTE: need to change the way hexagons are generated to accommodate
    #       instances where a user has an already pre-processed sfe object.
    sfe <- add.spotHex(sfe = sfe,
                       samples = samples,
                       sample_id = sample_id,
                       res = res,
                       flipped = flipped,
                       barcodes = barcodes)
  }

  ## Check and output either an msfe or an sfe object
  out <- .int_checkAndOutput(m_sfe = m_sfe, sfe = sfe, sample_id = sample_id)

  return(out)
}


