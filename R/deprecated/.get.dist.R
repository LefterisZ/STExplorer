#' @name .get.dist
#' 
#' @description A wrapper for base R dist function. Gives the ability to be run
#'              in a sapply function over a 3D array. The function is used by 
#'              the get.dist.array function.
#' @param data the input data, must be a 3D array.
#' @param focus The index for the location in focus.
#' @param method the distance measure to be used. This must be one of "euclidean", 
#'               "maximum", "manhattan", "canberra", "binary" or "minkowski". 
#'               Any unambiguous substring can be given. For more info, look at 
#'               the base R dist function.
#' @param p The power of the Minkowski distance. For more info, look at the base 
#'          R dist function.
#' 
#'

.get.dist <- function(focus, data, method = "euclidean", p = 2){
    d <- dist(data[,,focus], method, p) %>% 
        as.matrix()
    
    return(d)
}