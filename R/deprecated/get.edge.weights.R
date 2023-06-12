#' @name get.edge.weights
#' 
#' @description Build WEIGHTS vector for the graph edges. 
#' 
#' @param edges a two-column matrix containing the indexes of the two locations
#'              that make an edge (from --> to).
#' 
#' @param wdmat A distance matrix. It can be any distance matrix of size n x n
#'              (n = the total number of locations) but preferably this matrix
#'              should be generated with a distance decay function (i.e., 
#'              gw.weight) so that the distances between locations are weighted 
#'              and can be used as weights for the graph edges.
#'
#' @return a vector of weights matching the rows from the edges argument input.
#'
#' @export

get.edge.weights <- function(edges, wdmat){
    apply(edges, 1, function(x){
        wdmat[x[1],x[2]]
    })
}