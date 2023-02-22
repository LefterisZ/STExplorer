#' This script contains two pasting functions that are used internally by the 
#' .get.gwGraph function to generate the flags. 
#' 
#' These functions are wrapping the Sort and Match functions from the Rfast 
#' package. Rfast package implements basic R functions in C++, making them even
#' faster.
#' 
#' @param x the vector to paste together
#' 

## sort and collapse
.sort.N.paste <- function(x){
  paste0(sort(x), collapse = "_")
}

## collapse without sorting
.paste <- function(x){
  paste0(x, collapse = "_")
}
