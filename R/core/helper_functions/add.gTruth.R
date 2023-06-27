#' Add ground truth annotation to SpatialFeaturesExperiment object
#'
#' @name add.gTruth
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to add the ground truth annotation in the colData slot of a
#' SpatialFeaturesExperiment (SFE) object. Multiple annotations can be imported
#' using this approach, but note that if not all spots are annotated, NAs will
#' be imported. The SFE object has other slots for adding extra annotations that
#' may not cover the entire dataset.
#'
#' @param sfe The SpatialFeaturesExperiment object.
#' @param gtruth A dataframe containing the ground truth for the dataset.
#' It should have at least three columns: "Barcode" (spot barcodes matching
#' the colnames of the SFE object), "sample_id" (sample ID for each spot,
#' important for multiple samples), and the ground truth annotation column.
#'
#' @return The modified SpatialFeaturesExperiment object with the ground truth
#' annotation added to the colData slot.
#'
#' @details
#' This function adds the ground truth annotation to the colData slot of a
#' SpatialFeaturesExperiment object. It performs a merge between the colData of
#' the SFE object and the provided ground truth dataframe based on "Barcode" and
#' "sample_id". The merge includes all spots, and any missing annotations are
#' represented as NAs in the merged colData. The row names of the colData are
#' temporarily stored in a separate column before merging and then set back as
#' row names after the merge operation.
#'
#' @importFrom SpatialFeatureExperiment colData
#' @importFrom S4Vectors DataFrame
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso
#' \code{\link{add.index}}: A function to add an index to the colData slot of a
#' SpatialFeaturesExperiment object.
#' \code{\link{add.barcodes}}: A function to add barcodes to the colData slot of
#' a SpatialFeaturesExperiment object.
#'
#'
add.gTruth <- function(sfe, gtruth) {
  ## Add the row names in a column because after merging are being lost.
  colData(sfe)$rowNames <- colnames(sfe)
  ## Merge colData and ground truth
  colData(sfe) <- merge(colData(sfe), DataFrame(gtruth),
                        by = c("Barcode", "sample_id"), all = TRUE)
  ## Set rownames back
  rownames(colData(sfe)) <- colData(sfe)$rowNames
  colData(sfe)$rowNames <- NULL

  return(sfe)
}
