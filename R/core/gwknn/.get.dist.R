#' @name .get.dist
#' 
#' @description A wrapper for base R dist function. Gives the ability to run
#'              over a 3D array and to also return a 3D array.
#' @param data the input data (array, matrix, df), must be 2D.
#' @param focus 
#' 
#'

.get.dist <- function(data, focus, method = "euclidean", p = 2){
    d <- dist(data[,,focus], method, p) %>% 
        as.matrix()
    
    return(d)
}