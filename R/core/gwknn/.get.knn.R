#' @name get.knn
#' 
#' @description A function to find the k nearest neighbours
#' 
#' @param 
#' @param 
#' @param
#' 
#' @export

nearest_neighbors = function(x,obs, k, FUN, p = NULL){
    
    # Check the number of observations is the same
    if(ncol(x) != ncol(obs)){
        stop('Data must have the same number of variables')
    }
    
    # Calculate distance, considering p for Minkowski
    if(is.null(p)){
        dist = apply(x,1, FUN,obs)  
    }else{
        dist = apply(x,1, FUN,obs,p)
    }
    
    # Find closest neighbours
    distances = sort(dist)[1:k]
    neighbor_ind = which(dist %in% sort(dist)[1:k])
    
    if(length(neighbor_ind)!= k){
        warning(
            paste('Several variables with equal distance. Used k:',length(neighbor_ind))
        )
    }
    
    ret = list(neighbor_ind, distances)
    return(ret)
}