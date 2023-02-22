#' This script contains two sorting functions that are used internally by the 
#' .get.knn function. 
#' 
#' These functions are wrapping the Sort and Match functions from the Rfast 
#' package. Rfast package implements basic R functions in C++, making them even
#' faster.
#' 
#' @param x the vector to sort
#' @param part the number for partial sorting
#' @param .coln the colnames from the distance matrix that is provided to the 
#'              .get.knn function
#' 


.Sort.part <- function(x, part){
  Rfast::Sort(x, partial = part)[1:part]
}

.Match.Sort.part <- function(x, part, .coln) {
  .coln[match(Rfast::Sort(x, partial = part)[1:part], x)]
}
