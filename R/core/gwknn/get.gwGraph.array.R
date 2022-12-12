#' @name get.gwGraph.array
#' 
#' @description A function that uses a knn list object that contains indexes and
#'              distances of neighbours and transforms it into a 3D array of 
#'              graph matrices that contain 3 columns: "From", "To" and "W.Dist".
#' 
#' @param kList a list of lists of neighbour indexes and distances. (can be 
#'              created using the get.gwKNN.list function). More specifically, 
#'              it is a list that contains sub-lists of two matrices. Each 
#'              sub-list must have a matrix named "indexes" and a matrix named 
#'              "distances".
#'
#' @param focus.n numeric vector of indexes of locations to be used. If nothing
#'                is provided then all locations will be used.
#' 
#' @export

get.gwGraph.array <- function(kList, focus.n){
    
    # Check that the list is OK
    if(!is.list(kList)){
        stop("kList argument must be supplied with a list")
    }
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:length(kList) # indexes of locations to use (z-axis)
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of locations was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values only. Please make sure you provide location indexes only.")
    }
    
    # Get the sub-list names
    names <- names(kList)[focus.n]
    
    # Progress info
    message("Total number of graphs to generate: ", length(names))
    
    # Get the graphs
    temp <- sapply(seq_along(names), 
                   .get.gwGraph,
                   kList = kList,
                   names = names,
                   simplify = "array")
    
    return(temp)
}
