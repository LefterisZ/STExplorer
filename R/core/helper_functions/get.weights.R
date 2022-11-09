#' @name get.weights
#' 
#' @description 
#' 
#' @param 
#' 
#' @export

get.weights <- function(edges, wdmat){
    apply(edges, 1, function(x){
        wdmat[x[1],x[2]]
    })
}