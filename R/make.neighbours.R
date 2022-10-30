#' @name  make.neighbours
#' 
#' @description 
#' 
#' @param 
#' 
#' @export

make.neighbours <- function(){
    ## Create contiguity neighbours
    neighbours <- poly2nb(polygons, snap = 0)
    names(neighbours) = attr(neighbours, "region.id") # add names to the sub-lists
}
