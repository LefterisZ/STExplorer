#' @name get.gwGraph.table
#' 
#' @description a function to generate a k-NN graph from a gene expression 
#'              matrix. The gene expression is geographically weighted, 
#'              a statistical distance is calculated, k-Nearest Neighbours are
#'              identified and a graph is generated.
#' @param ..focus the index of the location in focus (needs to be character).
#' @param ..obs the gene expression matrix
#' @param ..wdmat the geographically weighted physical distance matrix
#' @param ..method the method to be used for the statistical distance. Default is 
#'               "euclidean".
#' @param ..k k for the k-NN. Keep in mind that a self neighbour is calculated and
#'          then dropped. So, if you want 7 neighbours you must set k = 8.
#' @param ..combine combine reciprocal edges. This means that A --> B is 
#'                 considered the same as B --> A. Default is FALSE.
#' @param ... passes arguments to the Rfast::Dist function.
#' @export 
#' 

get.gwGraph.table <- function(..focus, 
                              ..obs, 
                              ..wdmat, 
                              ..method, 
                              ..k,
                              ..combine,
                              ...) {
  
  # Get the network graph
  .get.gw.counts(focus = ..focus, obs = ..obs, wdmat = ..wdmat) %>% # weight gene expression
    .dist.wrap(gwCounts = ., method = ..method, ...) %>% # Calculate distances
    .get.knn(dists.X = ., k = ..k) %>% # Get nearest neighbours
    .get.gwGraph(kList = ., combine = ..combine) # Get the graphs
  
  
}