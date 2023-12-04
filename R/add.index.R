#' Add index to SpatialFeaturesExperiment object
#'
#' @name add.index
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to add the index in the colData slot of a
#' SpatialFeaturesExperiment (SFE) object. The index serves as a unique
#' identifier for every spot, even if there are multiple samples in the SFE
#' object.
#'
#' @keywords internal
#'
#' @param sfe The SpatialFeaturesExperiment object.
#'
#' @return The modified SpatialFeaturesExperiment object with the index added to
#' the colData slot.
#'
#' @details
#' This function adds an index column to the colData slot of a
#' SpatialFeaturesExperiment object. The index is generated as "spot_" followed
#' by a numeric sequence, starting from 1 and incrementing for each spot.
#' The index provides a unique identifier for every spot, ensuring each spot can
#' be uniquely identified even in the presence of multiple samples.
#'
#' @importFrom SpatialFeatureExperiment colData
#'
#' @seealso
#' \code{\link{add.barcodes}}: A function to add barcodes to the colData slot of
#' a SpatialFeaturesExperiment object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
add.index <- function(sfe) {
  ## Get the length of the SFE object
  len <- dim(colData(sfe))[1]
  ## Add the index
  colData(sfe)$index <- sprintf("spot_%d", seq(1:len))

  return(sfe)
}
