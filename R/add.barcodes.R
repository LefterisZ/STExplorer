#' Add barcodes in the colData
#'
#' @name add.barcodes
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to add the barcodes in the \code{colData} slot of a
#' \code{SpatialFeaturesExperiment} (SFE) object. In case of multiple samples
#' the \code{read10xVisiumSFE} function adds a suffix in the barcodes to make
#' them unique so that they can be used as col/rownames. This function takes
#' care of this by removing the suffix when transferring the barcodes to a
#' column.
#'
#' @keywords internal
#'
#' @param sfe The SFE object.
#'
#' @details
#' This function adds a new column named "Barcode" to the \code{colData} slot
#' of the given \code{SpatialFeaturesExperiment} (SFE) object. The function
#' moves the existing row names to this new column. If the barcodes in the
#' \code{colData} have a suffix in the format "-(0-9)-(0-9)", which is added
#' by the \code{read10xVisiumSFE} function to make them unique for multiple
#' samples, this function removes the suffix from the barcodes in the new
#' column. Additionally, the function extracts the capture area index from the
#' barcodes and adds it as a new column named "Capt_area" in the \code{colData}.
#'
#' @param sfe The SpatialFeaturesExperiment object.
#'
#' @returns The modified SpatialFeaturesExperiment object with added barcode
#' columns.
#' @importFrom SpatialFeatureExperiment colData
#' @importFrom SpatialFeatureExperiment colData<-
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso
#' \code{\link{read10xVisiumSFE}}
#'
add.barcodes <- function(sfe) {
  ## Move rownames to a column
  colData(sfe)$Barcode <- rownames(colData(sfe))
  ## Find which rows have a suffix
  bc_check <- grepl("-[0-9]-[0-9]$", colData(sfe)$Barcode)
  ## Remove the suffix from these rows in the new column
  colData(sfe)$Barcode[bc_check] <- gsub("-[0-9]$", "",
                                         colData(sfe)$Barcode[bc_check])
  ## Add the capture area index
  colData(sfe)$Capt_area <- gsub("^[A,T,G,C]*-", "",
                                 colData(sfe)$Barcode)

  return(sfe)
}
