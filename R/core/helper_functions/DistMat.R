#' Add Distance Matrix
#'
#' Calculate and store the distance matrix based on spatial coordinates using
#' the specified metric.
#' @name addDistMat
#'
#' @param sfe \code{SpatialFeatureExperiment} (SFE) object.
#' @param p Numeric scalar specifying the power parameter for the Minkowski
#' distance. See details below for more information.
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
addDistMat <- function(sfe, p, ...) {
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

#' Get Distance Matrix
#'
#' Retrieve the distance matrix based on the specified metric from the metadata
#' of a SpatialFeatureExperiment object.
#'
#' @param sfe \code{SpatialFeatureExperiment} (SFE) object.
#' @param dMetric Character string specifying the distance metric.
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
getDistMat <- function(sfe, dMetric) {
  dMat <- sfe@metadata$dMat[[dMetric]]
  return(dMat)
}
