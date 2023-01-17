#' @name get.edge.freq
#' 
#' @description A function to look inside a 3D gwGraph array and find the frequency 
#'              of every existing edge and calculate a weight. The returned data
#'              frame contains 5 columns: "from", "to", "flag", 
#'              "weights" and "count". Where "from" and "to" columns contain the 
#'              spot names that make an edge between them, "flag" is an 
#'              identifier to identify the edge pair (note here that A -> B and 
#'              B -> A are regarded as the same edge and have the same flag), 
#'              "weights" is the sum of the weighted distance between the 
#'              nodes of the edge (every time an edge is found amongst the 
#'              k-nearest neighbours, the weight it bears is added to the sum)  
#'              and finally, "count" is the number of times this flag is present 
#'              amongst the k-nearest neighbours. 
#' 
#' @param wGraph a 3D gwGraph array that contains all gwKNN graphs as created by
#'               the get.gwGraph.array() function.
#' @param focus.n Numeric or character. The indexes (numeric) or the names 
#'                (character) of the locations you want in focus. It can be a 
#'                vector of indexes from locations that are of interest. Default
#'                behaviour is for focus.n to be missing. This will result to  
#'                all locations being considered.
#' 
#' @return a data frame with dims = [all possible edges, 5]
#' 
#' @export


get.edge.freq <- function(wGraph, focus.n){
    
    # Get dimensions to make the data frame to store data
    nds <- attr(wGraph, "spot_names")[[1]] # get the node indexes
    
    # Make the data frame with all possible edges between the nodes
    message("Getting unique set of edges between nodes ...")
    edges <- .get.all.edges(nodes = nds) %>% 
        mutate(count = 0) %>% 
        mutate(weights = 0)
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:dim(wGraph)[3] # indexes of locations to use (z-axis)
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of location indexes was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
    } else if(is.vector(focus.n) & is.character(focus.n)) {
        message("A selection of location names was provided...")
        message("Locations with names: ", paste(focus.n, collpse = " "))
        focus.n <- match(focus.n, dimnames(wGraph)[[3]])
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values OR character values only. 
             Please make sure you provide location indexes only.")
    }
    
    # Iterate through the locations to be used and update the edge matrix
    for(n in focus.n){
        edges <- .update.edge.mtx(n = n, 
                                  wGraphs_arr = wGraph,
                                  edges_mtx = edges, 
                                  focus.n = focus.n)
    }
    
    return(edges)
}

