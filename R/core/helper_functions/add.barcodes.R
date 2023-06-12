#' @name add.barcodes
#' 
#' @description
#' A function to add the barcodes in the \code{colData} slot of a 
#' \code{SpatialFeaturesExperiment} (SFE) object. In case of multiple samples 
#' the \code{read10XVisiumSFE} function adds a suffix in the barcodes to make
#' them unique so that they can be used as col/rownames. This function takes 
#' care of this by removing the suffix when transferring the barcodes to a 
#' column. 
#' 
#' @param obj The SFE object.
#' 
#' @export

add.barcodes <- function(obj) {
  ## Move rownames to a column
  colData(obj)$Barcode <- rownames(colData(obj)) 
  ## Find which rows have a suffix
  bc_check <- grepl("-[0-9]-[0-9]$", colData(obj)$Barcode)
  ## Remove the suffix from these rows in the new column
  colData(obj)$Barcode[bc_check] <- gsub("-[0-9]$", "", 
                                         colData(obj)$Barcode[bc_check])
  ## Add the capture area index
  colData(obj)$Capt_area <- gsub("^[A,T,G,C]*-", "", 
                                 colData(obj)$Barcode)
  
  return(obj)
}
