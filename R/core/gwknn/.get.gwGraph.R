#' @name .get.gwGraph
#' 
#' @description A function to generate a 3 column graph data frame. Columns are
#'              "From", "To" and "W.Dist". It is meant to be used internally  
#'              from get.gwGraph.array().
#' 
#' @param focus.nm the name of the location in focus. A character string of 
#'                 length 1.
#' @param kList the list containing sub-lists of indexes and distances from kNNs.
#' 


.get.gwGraph <- function(focus.nm, kList){
    
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
    
    # Attach distances to graph df
    graph <- graph %>%
        mutate("W.Dist" = edgeW) %>% # add the weighted distances as a third column
        as.matrix()
    
    return(graph)
}
