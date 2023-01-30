#' @name get.gwKNN.list
#' 
#' @description A function to find the k-nearest neighbours of a location and 
#'              returns a list of neighbour indexes and distances. When selecting 
#'              a k, keep in mind that the self neighbour is also retrieved.
#' 
#' @param dists A 3D array containing the weighted distances. (can be generated
#'              using the get.dist.array function)
#' @param k The number of neighbours to select. Remember that the self-neighbour 
#'          is included. Therefore if you want 6 nearest neighbours, you should 
#'          select k = 7 since 1 of the 7 neighbours is the self-neighbour.
#' @param focus.n Numeric or character. The indexes (numeric) or the names 
#'                (character) of the locations you want in focus. It can be a 
#'                vector of indexes from locations that are of interest. Default
#'                behaviour is for focus.n to be missing. This will result to  
#'                all locations being considered.
#' @param strategy a future plan strategy for parallelisation. The defaults is 
#'                 multicore because it doesn't need to load the environment on 
#'                 each worker and thus is more suitable for large datasets. For
#'                 more info look at future package.
#' @param workers the number of cores to use for parallelisation.
#' 
#' @return a list with length = focus.n. Each list entry includes two matrices 
#'         one of the indexes (spot names) and one of the distances of the 
#'         k-neighbours.
#'         
#' @export

get.gwKNN.list <- function(dists, k, focus.n, strategy = "multicore", workers = 5){
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:dim(dists)[3] # indexes of locations to use (z-axis)
        names(focus.n) <- rownames(dists) # give spot names to the indexes
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of location indexes was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
        names(focus.n) <- dimnames(dists)[[1]][focus.n] # give names
    } else if(is.vector(focus.n) & is.character(focus.n)) {
        message("A selection of location names was provided...")
        message("Locations with names: ", paste(focus.n, collpse = " "))
        focus.n <- match(focus.n, dimnames(dists)[[3]]) # find the indexes
        names(focus.n) <- dimnames(dists)[[1]][focus.n] # give names
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values OR character values only. 
             Please make sure you provide location indexes only.")
    }
    
    # Find the neighbours without prallelisation
    # temp <- lapply(focus.n, 
    #                function(X, dists, k){
    #                    message("Getting k-NNs from distance matrix with index: ", X)
    #                    .get.knn(dists[,,X], k)
    #                },
    #                dists = dists,
    #                k = k)
    
    
    
    # Find the neighbours with future parallelisation
    # check if there are enough cores
    if(availableCores() <= workers) {
        message("You have less cores available than the number of workers you 
                asked for. Or you are asking to use all of your cores.")
        message("Number of cores available: ", availableCores())
        message("Number of workers asked: ", workers)
        message("Re-assinging the worker number to availableCores - 1")
        workers = availableCores() - 1
        message("Number of workers given: ", workers)
    }
    
    # set up strategy and progress handlers
    plan(strategy = strategy, workers = workers)
    handlers(global = TRUE)
    handlers("progress")
    
    message("Finding the Neighbours...")
    message("This might take some time depending on the number of spots you included.")
    
    knnProg_FUN <- function(focus.n, k){
        pr <- progressr::progressor(along = focus.n)
        future_lapply(focus.n,
                      .get.knn,
                      #dists = dists,
                      k = k,
                      future.globals = c("dists", "pr"))
    }
    
    temp <- knnProg_FUN(focus.n = focus.n,
                        k = k)
    
    # Set the indexes as names for the sub-lists
    names(temp) <- dimnames(dists)[[1]][focus.n]
    
    # Set an attirbute with all spot names 
    attr(temp, "spot_names") <- dimnames(dists)[1]
    
    message("step 3/6: DONE!!")
    
    return(temp)
}
