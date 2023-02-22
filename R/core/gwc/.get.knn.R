#' @name .get.knn
#' 
#' @description A function to find the k nearest neighbours.
#' 
#' @param dists.X distance matrix
#' @param k the number of neighbours
#' 

.get.knn <- function(dists.X, k) {

  ## Find closest neighbours
  ## a. get their distances
  distances <- Rfast::transpose(apply(dists.X, 1,
                                      .Sort.part,
                                      part = k))
  ## b. get their indexes
  coln <- colnames(dists.X)
  neighbour_idx <- t(apply(dists.X, 1,
                           .Match.Sort.part,
                           part = k,
                           .coln = coln))
  
  ## Prepare output
  out = list(neighbour_idx, distances)
  names(out) <- c("indexes", "distances")
  attr(out, "focus") <- attr(dists.X, "focus")
  
  return(out)
}

