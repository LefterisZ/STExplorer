#' @name get.gwKNN.list
#' 
#' @description A function to find the k-nearest neighbours of a location and 
#'              returns a list of neighbour indexes and distances. When selecting 
#'              a k, keep in mind that the self neighbour is also retrieved.
#' @param dists A 3D array containing the weighted distances. (can be generated
#'              using the get.dist.array function)
#' @param k The number of neighbours to select.
#' @param focus.n The indexes for the locations in focus. If not provided then 
#'                the function will use all locations.
#' 
#' @export

get.gwKNN.list <- function(dists, k, focus.n){
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:dim(dists)[3] # indexes of locations to use (z-axis)
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of locations was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values only. Please make sure you provide location indexes only.")
    }
    
    # Find the neighbours
    temp <- lapply(focus.n, 
                   function(X, dists, k){
                       message("Getting k-NNs from distance matrix with index: ", X)
                       .get.knn(dists[,,X], k)
                       },
                   dists = dists,
                   k = k)
    # Set the indexes as names for the sub-lists
    names(temp) <- dimnames(dists)[[1]][focus.n]
    
    return(temp)
}
