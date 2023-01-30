#' @name .get.knn
#' 
#' @description A function to find the k nearest neighbours. Is used by the 
#'              get.gwKNN.list function.
#' 
#' @param X index of the slice to be evaluated
#' @param k the number of neighbours
#' 


.get.knn <- function(X, k){
    
    # Get the slice from the array
    dists.X <- dists[,,X]
    
    # Find closest neighbours
    ## a. get their distances
    distances <- t(apply(dists.X, 1, 
                      function(x){
                          sort(x)[1:k]
                          }
                      ))
    ## b. get their indexes
    neighbour_idx <- t(apply(dists.X, 1,
                         function(x){
                             names(sort(x)[1:k])
                             }
                         ))
    ## c. add dimnames --> spot names
    # dimnames(distances)[[1]]
    # dimnames(neighbour_idx)[[1]]
    
    # Prepare output
    out = list(neighbour_idx, distances)
    names(out) <- c("indexes", "distances")
    
    pr()
    
    return(out)
}
