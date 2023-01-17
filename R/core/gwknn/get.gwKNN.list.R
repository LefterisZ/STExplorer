#' @name get.gwKNN.list
#' 
#' @description A function to find the k-nearest neighbours of a location and 
#'              returns a list of neighbour indexes and distances. When selecting 
#'              a k, keep in mind that the self neighbour is also retrieved.
#' 
#' @param dists A 3D array containing the weighted distances. (can be generated
#'              using the get.dist.array function)
#' @param k The number of neighbours to select. Remember that the self-neighbour 
#'          is included. Therefore if you want 6 nearest neighbours, you should 
#'          select k = 7 since 1 of the 7 neighbours is the self-neighbour.
#' @param focus.n Numeric or character. The indexes (numeric) or the names 
#'                (character) of the locations you want in focus. It can be a 
#'                vector of indexes from locations that are of interest. Default
#'                behaviour is for focus.n to be missing. This will result to  
#'                all locations being considered.
#' 
#' @return a list with length = focus.n. Each list entry includes two matrices 
#'         one of the indexes (spot names) and one of the distances of the 
#'         k-neighbours.
#'         
#' @export

get.gwKNN.list <- function(dists, k, focus.n){
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:dim(dists)[3] # indexes of locations to use (z-axis)
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of location indexes was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
    } else if(is.vector(focus.n) & is.character(focus.n)) {
        message("A selection of location names was provided...")
        message("Locations with names: ", paste(focus.n, collpse = " "))
        focus.n <- match(focus.n, dimnames(dists)[[3]])
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values OR character values only. 
             Please make sure you provide location indexes only.")
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
    
    # Set an attirbute with all spot names 
    attr(temp, "spot_names") <- dimnames(dists)[1]
    
    return(temp)
}
