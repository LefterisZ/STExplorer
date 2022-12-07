#' @name get.gwDist.array
#' 
#' @description A wrapper for base R dist function. Gives the ability to return
#'              a 3D array.
#' @param gwCounts the input data, must be a 3D array of gw-counts as generated
#'                 by the get.gwCount.array() function.
#' @param focus.n The indexes for the locations in focus. If not provided then 
#'                the function will use all locations.
#' @param method The distance measure to be used. This must be one of "euclidean", 
#'               "maximum", "manhattan", "canberra", "binary" or "minkowski". 
#'               Any unambiguous substring can be given. For more info, look at 
#'               the base R dist function. Default is euclidean.
#' @param p The power of the Minkowski distance. For more info, look at the base 
#'          R dist function. Default is 2 which is the euclidean distance.
#' 
#'
get.gwDist.array <- function(gwCounts, focus.n, method = "euclidean", p = 2){
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:dim(gwCounts)[3] # indexes of locations to use (z-axis)
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of locations was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values only. Please make sure you provide location indexes only.")
    }
    
    # Find the distances
    temp <- sapply(focus.n, 
                   function(X, method, p){
                       message("Calculating location: ", X)
                       dist(gwCounts[,,X], method, p) %>% as.matrix()
                       }, 
                   method = method, 
                   p = p, 
                   simplify = "array")
    
    ## add dimnames --> spot names as rows and columns
    dimnames(temp)[[1]] <- rownames(gwCounts)
    dimnames(temp)[[2]] <- rownames(gwCounts)
    dimnames(temp)[[3]] <- rownames(gwCounts)[focus.n]
    
    return(temp)
}
