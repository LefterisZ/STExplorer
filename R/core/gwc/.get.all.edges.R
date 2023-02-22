#' @name .get.all.edges
#' 
#' @description a helper function called internally by get.edge.freq(). Finds 
#'              all putative node combinations that can make an edge and stores
#'              only the unique ones. For example, pair 1-->2 and pair 2-->1 are
#'              concidered to be the same.
#' 
#' @param nodes a vector of node indexes.
#' 
#' 

.get.all.edges <- function(nodes){
    
    # Set a function for the apply()
    ## sort and collapse
    .paste <- function(x){
        paste0(x, collapse = "_")
    }
    
    # Get the unique edges
    edge.comb <- expand.grid(nodes, nodes, stringsAsFactors = FALSE) %>% # get all possible edges
        mutate(flag = apply(., 1, .paste)) # get a flag column
        #.[!duplicated(.$flag),] # remove those duplicates
    
    colnames(edge.comb) <- c("from", "to", "flag")
    
    return(edge.comb)
}

