#' @name .get.gwGraph
#' 
#' @description A function to generate a 4 column graph data frame. Columns are
#'              "From", "To", "W.Dist" and "count".
#' 
#' @param kList a list with indexes of neighbours and their distances from the 
#'              location in focus.
#' @param combine combine reciprocal edges. This means that A --> B is 
#'                 considered the same as B --> A. Default is FALSE.
#' 

.get.gwGraph <- function(kList, 
                         combine = FALSE) {
    
    ## Build graph df From -> To
    graph <- kList$indexes # select indexes table
    colnames(graph) <- c("from", paste0("nb", 1:(dim(graph)[2] - 1))) # add colnames
    graph <- graph %>%
        as.data.frame() %>%
        pivot_longer(-from, names_to = NULL, values_to = "to") # pivot long
    
    ## Get weighted distances
    edgeW <- kList$distances # select distances table
    edgeW <- edgeW %>%
        as.data.frame() %>%
        dplyr::select(-c("V1")) %>% # remove first column because is self-neighbour
        t() %>% # transpose df to a matrix
        as.vector() # make the matrix a single vector row-wise
    
    ## Attach a flag to each edge (the concatenation of the sorted node indexes)
    if (combine) {
      graph <- graph %>% 
        mutate(flag = apply(., 1, .sort.N.paste))
    } else {
      graph <- graph %>% 
        mutate(flag = apply(., 1, .paste))
    }
    
    ## Attach distances to graph df
    graph <- graph %>%
        mutate(wDist = edgeW,
               count = 1) %>% # add the weighted distances as a third column
        as.matrix()
    
    ## Set attribute focus
    attr(graph, "focus") <- attr(kList, "focus")
    
    return(graph)
}
