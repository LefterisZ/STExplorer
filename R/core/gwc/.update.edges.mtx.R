#' @name .update.edge.mtx
#' 
#' @description A function to take a matrix out of a list that contains kNN 
#'              graphs (namely matrices having columns with "from", "to", 
#'              "weights" and a "count" column with 1 count for each edge)
#' 
#' @param focus the index of a spot as character
#' @param wGraphs_list the graphs list
#' @param edges_mtx the matrix with all the edges
#' 
#' 

.update.edge.mtx <- function(focus, wGraphs_list, edges_mtx){
    
    ## Take out a matrix from the list of graphs
    wGraph_s <- as.data.frame(wGraphs_list[focus],
                              col.names = dimnames(wGraphs_list[focus])[[2]])
    
    ## Transform counts and wDist columns to numerics
    wGraph_s$count <- as.numeric(wGraph_s$count)
    wGraph_s$wDist <- as.numeric(wGraph_s$wDist)
    
    ## Make a matching vector for the flags in question
    x <- match(wGraph_s[,"flag"], edges_mtx$flag, nomatch = 0) 
    
    ## Subset and replace the count
    edges_mtx[x, "count"] <- edges_mtx[x, "count"] + wGraph_s$count
    
    ## Subset and replace the weights
    edges_mtx[x, "wDistSum"] <- edges_mtx[x, "wDistSum"] + wGraph_s$wDist
    
    
    return(edges_mtx)
}

