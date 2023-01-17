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
    
    # Summarise edge counts per flag
    c <- tapply(as.numeric(wGraph_s$count), wGraph_s$flag, sum) %>%
        data.frame(flag = names(.), count_y = ., row.names = NULL) 
    
    # Summarise edge weights per flag
    wD <- tapply(as.numeric(wGraph_s$wDist), wGraph_s$flag, sum) %>% 
        data.frame(flag = names(.), weights_y = ., row.names = NULL)
    
    # Merge count and weight sums by flag
    wDc <- merge(c, wD)
    
    # Keep only the "from", "to" and "flag" columns from the slice-matrix and 
    # remove duplicated flags because this creates problems with the next merge 
    # by either inserting duplicates again or NAs
    wGraph_s <- wGraph_s[!duplicated(wGraph_s$flag), c("from", "to", "flag")] 
    
    # Merge "from", "to" and "flag" with summed count and weights columns by flag
    wGraph_s <- merge(wGraph_s, wDc, all.x = FALSE, all.y = TRUE)
    
    # Make a matching vector for the flags in question
    x <- match(wGraph_s[,"flag"], edges_mtx$flag, nomatch = 0) 
    
    # Subset and replace the count
    edges_mtx[x, "count"] <- edges_mtx[x, "count"] + wGraph_s$count_y
    
    # Subset and replace the weights
    edges_mtx[x, "weights"] <- edges_mtx[x, "weights"] + wGraph_s$weights_y
    
    
    return(edges_mtx)
}
