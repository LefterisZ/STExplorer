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
        mutate(weights_x = 0)
    
    
    
    
}




wGraphh = as.data.frame(wGraph[,,1]) # take out a slice from the 3D array
colnames(wGraphh) <- dimnames(wGraph)[[2]] #give it colnames

c <- tapply(as.numeric(wGraphh$count), wGraphh$flag, sum) %>% # find counts sum per flag
    data.frame(flag = names(.), count_y = ., row.names = NULL) 

wD <- tapply(as.numeric(wGraphh$wDist), wGraphh$flag, sum) %>% #find weights sum per flag
    data.frame(flag = names(.), weights_y = ., row.names = NULL)

wDc <- merge(c, wD) # merge them together

wGraphh <- wGraphh[!duplicated(wGraphh$flag), c("from", "to", "flag")] # remove duplicate flags because this creates problems with the next merge with either inserting duplicates again or NAs

wGraphh <- merge(wGraphh, wDc, all.x = FALSE, all.y = TRUE)

# x <- c(edges$flag %in% wGraphh[,"flag"])
# edges$count <- edges$count + x

z <- sqldf("SELECT A.*, A.weights_x + B.weights_y AS weights FROM edges A LEFT JOIN wGraphh B ON A.flag = B.flag")




wGraph = graph.W
focus = 1
focus.n = 1:3
k = 7
x = c(1, 2)

rm(wGraph, wGraphh, wGraph_1, obs.W2, focus, focus.n, edges, edges.2, edges.3, x, nR, n, nodes, edge.comb, edges.count, nds, sort.N.paste,x,y,z,h,inpck,envir)
