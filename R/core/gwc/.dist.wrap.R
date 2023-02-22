#' @name .dist.wrap
#' 
#' @description a wrapper for the Rfast Dist function to be used for calculating
#'              statistical distances from a gene expression matrix that is 
#'              weighted from distance from a focus point.
#' 
#' @param gwCounts the gene expression matrix
#' @param method the method to calculate the distance (more info: Rfast::Dist)
#' @param p the Minkowski power (more info: Rfast::Dist)
#' @param ... more arguments to be passed to Rfast::Dist (more info: Rfast::Dist)
#' 
#' @return a n x n matrix (n = number of locations)
#' 

.dist.wrap <- function(gwCounts, method, ...){
  
  ## Calculate distances
  out <- Rfast::Dist(gwCounts, method = method, ...)
  
  ## Add dimnames --> spot indexes as row and column names
  dimnames(out)[[1]] <- rownames(gwCounts)
  dimnames(out)[[2]] <- rownames(gwCounts)
  attr(out, "focus") <- attr(gwCounts, "focus")
  
  return(out)
}

