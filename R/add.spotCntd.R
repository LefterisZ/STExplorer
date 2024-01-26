#' Generate centroids from spot coordinates in a SpatialFeaturesExperiment
#' object
#'
#' @name add.spotCntd
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to generate centroids from the spot coordinates. It takes the
#' coordinates out of the \code{spatialCoords} slot of a
#' SpatialFeaturesExperiment (SFE) object and stores it inside the colGeometries
#' slot.
#'
#' @keywords internal
#'
#' @param sfe The SpatialFeaturesExperiment object.
#' @param sample_id A character string specifying unique sample identifiers, one
#' or each directory specified via \code{samples} when you loaded the SFE object
#' using the \code{read10xVisiumSFE} function.
#'
#' @details
#' This function extracts the spot coordinates from the \code{spatialCoords}
#' slot of the SFE object, generates centroids from the coordinates, and stores
#' the resulting centroids in the \code{colGeometries} slot of the SFE object.
#' The centroid coordinates are represented as an \code{sf} object, with each
#' centroid corresponding to a spot in the SFE object. The row names of the
#' \code{cntds} object are set to match the column names of the SFE object for
#' easy reference and linking between the two.
#'
#' @importFrom SpatialFeatureExperiment spatialCoords colGeometry
#' @importFrom SpatialFeatureExperiment colGeometry<-
#' @importFrom sf st_as_sf
#' @importFrom magrittr %>%
#'
#' @seealso
#' \code{\link{add.index}}: A function to add an index in the colData slot of a
#' SpatialFeaturesExperiment object.
#' \code{\link{add.barcodes}}: A function to add barcodes in the colData slot of
#' a SpatialFeaturesExperiment object.
#' \code{\link{add.gTruth}}: A function to add ground truth annotation in the
#' colData slot of a SpatialFeaturesExperiment object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @returns The modified SpatialFeaturesExperiment object with the centroid
#' coordinates stored in the colGeometries slot.
#'
add.spotCntd <- function(sfe, sample_id) {
  ## Fetch coordinates
  coords <- colnames(spatialCoords(sfe))
  ## Get centroids
  cntds <- spatialCoords(sfe) %>%
    as.data.frame() %>%
    st_as_sf(coords = coords,
             remove = TRUE)
  ## Add rownames
  rownames(cntds) <- colnames(sfe)
  ## Add inside colGeometries slot
  for (sID in sample_id) {
    selection <- colData(sfe)$sample_id %in% sID
    colGeometry(sfe, sample_id = sID, "spotCntd") <- cntds[selection, ]
  }

  return(sfe)
}
