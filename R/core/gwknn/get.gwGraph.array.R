#' @name get.gwGraph.array
#' 
#' @description A function that uses a knn list object that contains indexes and
#'              distances of neighbours and transforms it into a 3D array of 
#'              graph matrices that contain 3 columns: "From", "To" and "W.Dist".
#' 
#' @param kList a list of lists of neighbour indexes and distances. (can be 
#'              created using the get.gwKNN.list function). More specifically, 
#'              it is a list that contains sub-lists of two matrices. Each 
#'              sub-list must have a matrix named "indexes" and a matrix named 
#'              "distances".
#' 
#' @export

get.gwGraph.array <- function(kList){
    
    # Check that the list is OK
    if(!is.list(kList)){
        stop("kList argument must be supplied with a list")
    }
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:nrow(obs) # indexes of locations to use (z-axis)
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of locations was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values only. Please make sure you provide location indexes only.")
    }
    
    # Get the sub-list names
    names <- names(kList)
    
    # Get the graphs
    temp <- sapply(names, 
                   function(X, kList){
                       
                   },
                   kList = kList)
    
    return(temp)
}


graph <- kList[[1]]$indexes  %>% t()
colnames(graph) <- c("from", paste0("nb", 1:(dim(graph)[2]-1)))
graph <- graph %>%
    as.data.frame() %>%
    pivot_longer(-from, names_to = NULL, values_to = "to")

edgeW <- apply(graph, 1, 
               function(x){
                   dist.Mat.w1[x[1],x[2]]
                   }
               )

dists = dist.W[,,1]
kList = knn.W
focus.n <- 1:3
dists = dist.W[,,1]
k = 7
rm(dists, k, out, outer,kList,focus.n)
