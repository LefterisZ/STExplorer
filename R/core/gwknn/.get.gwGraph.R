#' @name .get.gwGraph
#' 
#' @description A function to generate a 3 column graph data frame. Columns are
#'              "From", "To" and "W.Dist". It is meant to be used internally  
#'              from get.gwGraph.array().
#' 
#' @param n a number
#' @param kList the list containing sub-lists of indexes and distances from kNNs.
#' 
#' @param names a vector with the names of the locations in focus.
#' 
#' 
# n = 1
# kList = knn.W
# names = names(knn.W)
# k = 7
# dists =dist.W[,,1]
# wGraph = graph.W
# rm(wGraph, nodes, n, kList, names, focus.nm,graph,edgeW,dists,distances, neighbour_idx,focus.n,temp,out)

.get.gwGraph <- function(n, kList, names){
    
    # Set a function for the apply()
    ## sort and collapse
    sort.N.paste <- function(x){
        paste0(sort(x), collapse = "_")
    }
    
    # Print a message to show progress
    message("Generating graph ", n, "/", length(names))
    
    # Get the name of the location
    focus.nm <- names[n]
    
    # Build graph df From -> To
    graph <- kList[[focus.nm]]$indexes # select indexes table
    colnames(graph) <- c("from", paste0("nb", 1:(dim(graph)[2]-1))) # add colnames
    graph <- graph %>%
        as.data.frame() %>%
        pivot_longer(-from, names_to = NULL, values_to = "to") # pivot long
    
    # Get weighted distances
    edgeW <- kList[[focus.nm]]$distances # select distances table
    edgeW <- edgeW %>%
        as.data.frame() %>%
        dplyr::select(-c("V1")) %>% # remove first column because is self-neighbour
        t() %>% # transpose df to a matrix
        as.vector() # make the matrix a single vector row-wise
    
    # Attach a flag to each edge (the concatenation of the sorted node indexes)
    graph <- graph %>% 
        # transform(from = as.integer(from), 
        #           to = as.integer(to)) %>% # transform to integers to sort arithmetically
        mutate(flag = apply(., 1, sort.N.paste))
    
    # Attach distances to graph df
    graph <- graph %>%
        mutate(wDist = edgeW,
               count = 1) %>% # add the weighted distances as a third column
        as.matrix()
    
    return(graph)
}
