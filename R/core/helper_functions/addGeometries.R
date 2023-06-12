#' @name addGeometries
#' 
#' @description
#' A function to add centroid and hexagon geometries in the \code{colGeometries}
#' slot of an SpatialFeatureExperiment (SFE) object.
#' 
#' @param obj The SpatialFeaturesExperiment object.
#' 
#' @param samples 	a character vector specifying one or more directories, each 
#' corresponding to a 10x Genomics Visium sample (see Details); if provided, 
#' names will be used as sample identifiers.
#' 
#' @param sample_id character string specifying unique sample identifiers, one 
#' for each directory specified via samples; ignored if 
#' !is.null(names(samples)).
#' 
#' @param res the desired resolution. Can take one of "lowres", "hires", 
#' "fullres".
#' 
#' @export

addGeometries <- function(obj,
                          samples,
                          sample_id,
                          res = c("lowres", "hires", "fullres")) {
  
  res <- match.arg(res)
  
  ## Add Centroids
  obj <- add.spotCntd(obj)
  
  ## Get/ calculate spot diameter
  obj <- spot_diameter(obj = obj, 
                       samples = samples, 
                       sample_id = sample_id, 
                       res = res)
  
  ## Add Hexagons
  obj <- add.spotHex(obj = obj, 
                     samples = samples, 
                     sample_id = sample_id, 
                     res = res)
  
  return(obj)
}
