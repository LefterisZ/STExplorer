#' @name .get.gwGraph
#' 
#' @description A function to generate a 3 column graph data frame. Columns are
#'              "From", "To" and "W.Dist". It is meant to be used internally  
#'              from get.gwGraph.array().
#' 
#' @param X a number. The location in focus.
#' 
#' 

.get.gwGraph <- function(X){
    
    # Set a function for the apply()
    ## sort and collapse
    .paste <- function(x){
        paste0(x, collapse = "_")
    }
    
    # Build graph df From -> To
    graph <- kList[[X]]$indexes # select indexes table
    colnames(graph) <- c("from", paste0("nb", 1:(dim(graph)[2]-1))) # add colnames
    graph <- graph %>%
        as.data.frame() %>%
        pivot_longer(-from, names_to = NULL, values_to = "to") # pivot long
    
    # Get weighted distances
    edgeW <- kList[[X]]$distances # select distances table
    edgeW <- edgeW %>%
        as.data.frame() %>%
        dplyr::select(-c("V1")) %>% # remove first column because is self-neighbour
        t() %>% # transpose df to a matrix
        as.vector() # make the matrix a single vector row-wise
    
    # Attach a flag to each edge (the concatenation of the sorted node indexes)
    graph <- graph %>% 
        mutate(flag = apply(., 1, .paste))
    
    # Attach distances to graph df
    graph <- graph %>%
        mutate(wDist = edgeW,
               count = 1) %>% # add the weighted distances as a third column
        as.matrix()
    
    pr()
    
    return(graph)
}
