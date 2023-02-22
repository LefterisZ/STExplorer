#' @name get.edge.freq
#' 
#' @description A function to look inside a list of graphs and find the frequency 
#'              of every existing edge and calculate a weight. The returned data
#'              frame contains 6 columns: "from", "to", "flag", "wDistSum", 
#'              "obsFreq" and "count". Where "from" and "to" columns contain the 
#'              spot names that make an edge between them, "flag" is an 
#'              identifier to identify the edge pair (note here that A -> B and 
#'              B -> A are NOT regarded as the same edge and have a different 
#'              flag), "wDistSum" is the sum of the weighted distance between  
#'              the nodes of the edge (every time an edge is found amongst the 
#'              k-nearest neighbours, the weight it bears is added to its sum),
#'              "obsFreq" is the observed frequency of this edge in the data set 
#'              and finally, "count" is the number of times this flag is present 
#'              amongst the k-nearest neighbours.
#' 
#' @param wGraph a list of graphs that contains all gwKNN graphs as created by
#'               the get.gwGraphs() function.
#' @param remove_0 TRUE or FALSE. Remove entries in the the final matrix with a 
#'                 count sum of ZERO. Defaults to TRUE and better be left like 
#'                 that.
#' @param verbose show progress bar? Default is TRUE.
#' 
#' @return a data frame with dims = [all possible edges, 6]
#' 
#' @export


get.edge.freq <- function(wGraph, remove_0 = TRUE, verbose = TRUE){
    
    ## Get dimensions to make the data frame to store data
    nds <- names(wGraph) # get the node indexes
    
    ## Make the data frame with all possible edges between the nodes
    message("Getting unique set of edges between nodes ...")
    edges <- .get.all.edges(nodes = nds) %>% 
        mutate(count = 0) %>% 
        mutate(wDistSum = 0)
    
    ## Iterate through the locations to be used and update the edge matrix
    if (verbose) {
        message("Started to update the matrix. Progress notifications are set to every 100 iterations")
        pr <- progressr::progressor(along = nds)
        for (focus in nds) {
            edges <- .update.edge.mtx(focus = focus, 
                                      wGraphs_list = wGraph,
                                      edges_mtx = edges)
            if (as.numeric(focus) %% 10 == 0) {
                pr(message = sprintf("Updating edge matrix... %s / %d", 
                                     focus, length(nds)))
            }
        }
    } else {
        for (focus in nds) {
            edges <- .update.edge.mtx(focus = focus, 
                                      wGraphs_list = wGraph,
                                      edges_mtx = edges)
        }
    }
    
    
    ## Remove redundant edges that have a count of zero
    message("removing redundants ...")
    if (remove_0) {
        edges <- filter(edges, count > 0)
    }
    
    # Calculate the observed frequency of each edge
    edges <- mutate(edges, obsFreq = count/max(count))
    
    return(edges)
}

