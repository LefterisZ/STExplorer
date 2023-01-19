#' @name .update.edge.mtx
#' 
#' @description A function to take a slice out of a 3D array that contains kNN 
#'              graphs (namely matrices having columns with "from", "to", 
#'              "weights" and a "count" column with 1 count for each )
#' 
#' @param n a number
#' @param wGraphs_arr the graphs 3D array
#' @param edges_mtx the matrix with all the edges
#' @param focus.n a vector with the index of the locations in focus.
#' 
#' 

.update.edge.mtx <- function(n, wGraphs_arr, edges_mtx, focus.n){
    
    # Print a message to show progress
    message("Updating graph ", n, "/", length(focus.n))
    
    # Get the index of the location
    focus <- focus.n[n]
    
    # Take out a slice from the 3D array and give it colnames
    wGraph_s <- as.data.frame(wGraphs_arr[,,focus]) 
    colnames(wGraph_s) <- dimnames(wGraphs_arr)[[2]]
    
    # Transform counts and wDist columns to numerics
    wGraph_s$count <- as.numeric(wGraph_s$count)
    wGraph_s$wDist <- as.numeric(wGraph_s$wDist)
    
    # Make a matching vector for the flags in question
    x <- match(wGraph_s[,"flag"], edges_mtx$flag, nomatch = 0) 
    
    # Subset and replace the count
    edges_mtx[x, "count"] <- edges_mtx[x, "count"] + wGraph_s$count
    
    # Subset and replace the weights
    edges_mtx[x, "wDistSum"] <- edges_mtx[x, "wDistSum"] + wGraph_s$wDist
    
    
    return(edges_mtx)
}
