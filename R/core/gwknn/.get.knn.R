#' @name .get.knn
#' 
#' @description A function to find the k nearest neighbours. Is used by the 
#'              get.gwKNN.list function.
#' 
#' @param dists a matrix of distances
#' @param k the number of neighbours
#' 
#' @export

.get.knn <- function(dists, k){
    
    # Find closest neighbours
    # a. get their distances
    distances <- t(apply(dists, 1, 
                      function(x){
                          sort(x)[1:k]
                          }
                      ))
    # b. get their indexes
    neighbor_ind <- apply(dists, 1,
                         function(x){
                             which(x %in% sort(x)[1:k])
                             }
                         )
    
    
    # if(length(neighbor_ind)!= k){
    #     warning(
    #         paste('Several variables with equal distance. Used k:',length(neighbor_ind))
    #     )
    # }
    
    out = list(neighbor_ind, distances)
    names(out) <- c("indexes", "distances")
    
    
    
    return(out)
}
