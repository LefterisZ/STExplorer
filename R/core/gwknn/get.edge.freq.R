#' @name get.edge.freq
#' 
#' @description A function to look inside a 3D gwGraph aray and find the frequency 
#'              of every existing edge and calculate a weight.
#' 
#' @param wGraph a 3D gwGraph array that contains all gwKNN graphs as created by
#'               the get.gwGraph.array() function.
#' @param focus.n the index of the locations in focus. If left empty then all 
#'                locations are used.
#' @param k the k used for getting the k-NNs.
#' 
#' @export


get.edge.freq <- function(wGraph, focus.n, k){
    
    # Get dimensions to make the data frame to store data
    ## k-1 because we don't account for self-neighbours
    n <- dim(wGraph)[1]/(k-1) # get the number of nodes
    nds <- seq_len(n) # get the node indexes
    
    # Make the data frame
    message("Getting unique set of edges between nodes ...")
    edges <- .get.all.edges(nodes = nds) %>% 
        mutate(count = 0) %>% 
        mutate(weights = 0)
    
    
    
    
}


# edges <- wGraph[,,focus] %>%
#     as.data.frame() %>% 
#     filter(from == focus) %>% 
#     mutate("W.Dist" = as.numeric(W.Dist))
# 
# edges.2 <- as.matrix(wGraph[,,])

# edges.3 <- apply(wGraph, 2, c) # stack all matrices on top of each other
# # edges.3 <- apply(wGraph, 2, identity) # stack them maybe a bit faster
# edges.3 <- edges.3 %>% 
#     as.data.frame() %>% 
#     group_by(from, to)
# 
# edges.count <- count(edges.3) # count how many members each group has
# 
# edges.meanWdist <- edges.3 %>% 
#     mutate(meanWdist = mean(W.Dist))
# edge.comb <- expand.grid(nodes, nodes, stringsAsFactors = FALSE) %>% # get all possible edges
#     unite("pair", Var1, Var2, sep = ",", remove = TRUE) %>% 
#     mutate(pair = strsplit(.data$pair,","), .keep = "none")
# 
# mutate(temp = apply(.data, 1, 
#                     function(x){
#                         paste0(sort(x), collapse = "")
#                     }))
# apply(., 1, 
#       function(x){
#           paste0(sort(x), collapse = "")
#       }
# )

wGraph = graph.W
wGraph.1 = wGraph[,,1]
focus = 1
focus.n = 1:3
k = 7
x = c(1, 2)

rm(wGraph, obs.W2, focus, focus.n, edges, edges.2, edges.3, x, nR, n, nodes, edge.comb, edges.count, nds, sort.N.paste)
