#' @name add.index
#' 
#' @description
#' A function to add the index in the \code{colData} slot of a 
#' \code{SpatialFeaturesExperiment} (SFE) object. The index is a unique 
#' identifier for every spot even if there are multiple samples in the 
#' SFE object. 
#' 
#' @param obj The SFE object.
#' 
#' @export

add.index <- function(obj) {
  ## Get the length of the SFE object
  len <- dim(colData(obj))[1]
  ## Add the index
  colData(obj)$index <- sprintf("spot_%d", seq(1:len))
  
  return(obj)
}
