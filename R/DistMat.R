#' Add Distance Matrix
#'
#' Calculate and store the distance matrix based on spatial coordinates using
#' the specified metric.
#' @name addDistMat
#'
#' @param msfe \code{MetaSpatialFeatureExperiment} (Meta-SFE) object.
#' @param p Numeric scalar specifying the power parameter for the Minkowski
#' distance. See details below for more information.
#' @param sample_id character string or TRUE specifying sample/image
#' identifier(s); here, TRUE is equivalent to all samples/images.
#' @param ... Additional arguments to be passed to \code{gw.dist} function.
#'
#' @importFrom GWmodel gw.dist
#'
#' @return An updated SpatialFeatureExperiment object with the distance matrix
#' stored in the metadata.
#'
#' @details
#' This function calculates the distance matrix based on the spatial coordinates
#' of the provided SpatialFeatureExperiment object. The distance matrix is
#' computed using the specified metric, where p = 2 corresponds to Euclidean
#' distance, p = 1 corresponds to Manhattan distance, and p > 2 corresponds to
#' Minkowski distance with the specified power. The resulting distance matrix is
#' stored in the metadata of the SpatialFeatureExperiment object.
#'
#' @seealso [getDistMat()] for a function to get the distance matrix from within
#' an SFE object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Load sfe object
#' data(sfe)
#'
#' # Add distance matrix to sfe
#' sfe <- addDistMat(sfe, p = 2)
#'
#'
#' @export
addDistMat <- function(msfe, p, sample_id = TRUE, ...) {
  ## check arguments
  # stopifnot(is(msfe, "SpatialFeatureExperiment"))

  ## Select samples
  ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)

  ## Generate the graphs
  msfe_int <- lapply(msfe[ids], .int_addDistMat, p = p, ...)

  ## If specific samples where modified replace in the metaSFE list
  if (is.character(sample_id)) {
    msfe[names(msfe_int)] <- msfe_int
  } else {
    msfe <- msfe_int
  }

  return(msfe)
}

#' Get Distance Matrix
#'
#' Retrieve the distance matrix based on the specified metric from the metadata
#' of a SpatialFeatureExperiment object.
#'
#' @param msfe \code{MetaSpatialFeatureExperiment} (Meta-SFE) object.
#' @param dMetric Character string specifying the distance metric.#'
#' @param sample_id character string, or TRUE specifying sample/image
#' identifier(s); here, TRUE is equivalent to all samples/images.
#'
#' @return Distance matrix corresponding to the specified metric.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Load sfe object
#' data(sfe)
#'
#' # Add distance matrix to sfe
#' sfe <- addDistMat(sfe, p = 2)
#'
#' # Fetch distance matrix from sfe
#' dMat <- getDistMat(sfe, dMetric = "euclidean")
#'
#' # Check distance matrix
#' dMat[1:4, 1:4]
#'
#' @export
getDistMat <- function(msfe, dMetric, sample_id = TRUE) {
  ## Check sfe or msfe
  if (is(msfe, "SpatialFeatureExperiment")) {
    sfe <- TRUE
  } else if (is(msfe, "MetaSpatialFeatureExperiment")) {
    sfe <- FALSE
  }

  if (sfe) {
    dMat_int <- .int_getDistMat(msfe, dMetric = dMetric)
  } else if (!sfe) {
    ## Select samples
    ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)
    dMat_int <- lapply(msfe@sfe_data[ids], .int_getDistMat, dMetric = dMetric)
  }

  return(dMat_int)
}


#' Internal Function: Add Distance Matrix to SpatialFeatureExperiment
#'
#' This function calculates the distance matrix for spatial coordinates and
#' adds it to the metadata of a SpatialFeatureExperiment object.
#'
#' @param sfe A SpatialFeatureExperiment object containing spatial coordinates.
#' @param p Numeric, the power parameter for the distance metric (default is 2
#' for Euclidean distance).
#' @param ... Additional arguments to be passed to the gw.dist function.
#'
#' @return A modified SpatialFeatureExperiment object with the distance matrix
#' added to its metadata.
#'
#' @details The function calculates the distance matrix for the spatial
#' coordinates in the SpatialFeatureExperiment object using the specified
#' distance metric (Euclidean, Manhattan, or Minkowski). The resulting distance
#' matrix is stored in the metadata of the SpatialFeatureExperiment object.
#'
#' @seealso \code{\link{gw.dist}}, \code{\link[spatialCoords]{spatialCoords}},
#' \code{\link[colData]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords distance matrix, spatial coordinates, SpatialFeatureExperiment
#'
#' @rdname dot-int_addDistMat
#'
#' @importFrom GWmodel gw.dist
#'
.int_addDistMat <- function(sfe, p, ...) {
  dMat <- gw.dist(spatialCoords(sfe), p = p, ...)
  dimnames(dMat) <- list(NULL, colData(sfe)$Barcode)

  # Determine the name of the distance metric
  if (p == 2) {
    dMetric <- "euclidean"
  } else if (p == 1) {
    dMetric <- "manhattan"
  } else if (p > 2) {
    dMetric <- paste0("minkowski_", p)
  }

  # Store the distance matrix in the metadata of the SFE object
  sfe@metadata$dMat[[dMetric]] <- dMat

  return(sfe)
}

#' Internal Function: Get Distance Matrix from SpatialFeatureExperiment
#'
#' This function retrieves the distance matrix stored in the metadata of a
#' SpatialFeatureExperiment object.
#'
#' @param sfe A SpatialFeatureExperiment object containing a distance matrix in
#' its metadata.
#' @param dMetric The name of the distance metric used for the stored matrix
#' (Euclidean, Manhattan, or Minkowski).
#'
#' @return The distance matrix retrieved from the metadata of the
#' SpatialFeatureExperiment object.
#'
#' @details The function retrieves the distance matrix stored in the metadata
#' of the SpatialFeatureExperiment object. Users need to specify the distance
#' metric (Euclidean, Manhattan, or Minkowski) used for the stored matrix.
#'
#' @seealso \code{\link[colData]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords distance matrix, metadata, SpatialFeatureExperiment
#'
#' @rdname dot-int_getDistMat
#'
.int_getDistMat <- function(sfe, dMetric) {
  dMat <- sfe@metadata$dMat[[dMetric]]

  return(dMat)
}
