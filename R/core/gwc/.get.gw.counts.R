#' @name .get.gw.counts
#' 
#' @description A function that weights gene expression data (counts) based on a 
#'              focus location and on a distance matrix created using a preferable 
#'              kernel.
#' 
#' @param obs A matrix containing the observation data, usually gene expression 
#'            data where columns are the genes and rows are the locations.
#' @param wdmat An n x n matrix (n = number of locations) of weighted distances 
#'              between all locations.
#' @param focus The row name of the location in focus
#' 
#' @return a weighted gene expression data matrix
#' 
#'


.get.gw.counts <- function(obs, focus, wdmat){
    
    # Select weights for focus point
    w.focus <- wdmat[focus,]
    sumW <- sum(w.focus)                       ## get the weights sum for focus
    sqrt.W <- sqrt(w.focus)                    ## get sqrt of weights for focus
    
    # Apply weights
    w.obs.1 <- sweep(obs, 1, w.focus, "*") ## weight observations (1)
    colS <- colSums(w.obs.1)                   ## get the weighted colsums per gene
    w.obs.2 <- sweep(obs, 2, colS/sumW)        ## weight observations (2)
    w.obs <- sweep(w.obs.2, 1, sqrt.W, "*")## weight observations (final)
    
    w.obs <- as.matrix(w.obs)

    # add the row name of the location in focus as an attribute
    attr(w.obs, "focus") <- focus
    
    return(w.obs)
}
