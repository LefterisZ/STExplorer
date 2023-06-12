#' @name add.gTruth
#' 
#' @description
#' A function to add the ground truth annotation in the \code{colData} slot of a 
#' \code{SpatialFeaturesExperiment} (SFE) object. Of course multiple annotations
#' can be imported with this approach but keep in mind that if not all spots are
#' annotated then NAs will be imported. The \code{SFE} object has other slots 
#' that are for adding extra annotations that are not always covering the whole
#' dataset.
#' 
#' @param obj The SFE object.
#' @param gtruth A dataframe that contains the ground truth for your dataset. 
#' It needs to have at least 3 columns. One column named "Barcode" with the 
#' spot Barcodes (these need to match the colnames of the SFE object), one 
#' column named "sample_id" which includes the sample ID for each spot (this is
#' important when you import multiple samples) and one column with the ground 
#' truth annotation.
#' 
#' @export

add.gTruth <- function(obj, gtruth) {
  ## Add the row names in a column because after merging are being lost.
  colData(obj)$rowNames <- colnames(obj)
  ## Merge colData and ground truth
  colData(obj) <- merge(colData(obj), DataFrame(gtruth), 
                        by = c("Barcode", "sample_id"), all = TRUE)
  ## Set rownames back
  rownames(colData(obj)) <- colData(obj)$rowNames
  colData(obj)$rowNames <- NULL
  
  return(obj)
}
