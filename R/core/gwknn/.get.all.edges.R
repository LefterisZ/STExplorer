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
    sort.N.paste <- function(x){
        paste0(sort(x), collapse = "")
    }
    
    # Get the unique edges
    edge.comb <- expand.grid(nodes, nodes, stringsAsFactors = FALSE) %>% # get all possible edges
        mutate(temp = apply(., 1, sort.N.paste)) %>% # get a temp column that will have duplicates
        .[!duplicated(.$temp),] # remove those duplicates
    
    return(edge.comb)
}

