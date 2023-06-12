#' @name add.spotCntd
#' 
#' @description
#' A function to generate centroids from the spot coordinates. It takes the 
#' coordinates out of the \code{spatialCoords} slot of a 
#' \code{SpatialFeaturesExperiment} (SFE) object and stores it inside the 
#' \code{colGeometries} slot.
#' 
#' @param obj The SFE object.
#' 
#' @export

add.spotCntd <- function(obj) {
  ## Fetch coordinates
  coords <- colnames(spatialCoords(obj))
  ## Get centroids
  cntds <- spatialCoords(obj) %>% 
    as.data.frame() %>% 
    st_as_sf(coords = coords, 
             remove = TRUE)
  ## Add rownames
  rownames(cntds) <- colnames(obj)
  ## Add inside colGeometries slot
  colGeometry(obj, "spotCntd") <- cntds
  
  return(obj)
}
