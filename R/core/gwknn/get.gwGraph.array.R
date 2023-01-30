#' @name get.gwGraph.array
#' 
#' @description A function that uses a knn list object that contains indexes and
#'              distances of neighbours and transforms it into a 3D array of 
#'              graph matrices that contain 5 columns: "from", "to", "flag", 
#'              "wDist" and "count". Where "from" and "to" columns contain the 
#'              spot names that make an edge between them, "flag" is an 
#'              identifier to identify the edge pair (note here that A -> B and 
#'              B -> A are regarded as the same edge and have the same flag), 
#'              "wDist" is the weighted distance between the nodes of the edge 
#'              that is also used as an edge weight) and finally, "count" is the 
#'              number of times this flag is present amongst the k-nearest 
#'              neighbours.
#' 
#' @param kList a list of lists of neighbour indexes and distances. (can be 
#'              created using the get.gwKNN.list function). More specifically, 
#'              it is a list that contains sub-lists of two matrices. Each 
#'              sub-list must have a matrix named "indexes" and a matrix named 
#'              "distances".
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
#' 
#' @return a 3D array with dims = [graph edges, 5, focus.n]
#' 
#' @export

get.gwGraph.array <- function(kList, focus.n, strategy = "multicore", workers = 5){
    
    # Check that the list is OK
    if(!is.list(kList)){
        stop("kList argument must be supplied with a list")
    }
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:length(kList) # indexes of locations to use (z-axis)
        names <- names(kList)[focus.n] # get the sub-list names
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of location indexes was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
        names(focus.n) <- names(kList)[focus.n] # give names
    } else if(is.vector(focus.n) & is.character(focus.n)) {
        message("A selection of location names was provided...")
        message("Locations with names: ", paste(focus.n, collpse = " "))
        focus.n <- match(focus.n, names(kList)) # find the indexes
        names(focus.n) <- names(kList)[focus.n] # give names
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values OR character values only. 
             Please make sure you provide location indexes only.")
    }
    
    # Progress info
    message("Total number of graphs to generate: ", length(names))
    
    # Get the graphs without prallelisation
    # temp <- sapply(seq_along(names), 
    #                .get.gwGraph,
    #                kList = kList,
    #                names = names,
    #                simplify = "array")
    
    # Get the graphs with future parallelisation
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
    
    message("Generating the graphs...")
    
    graphProg_FUN <- function(focus.n){
        pr <- progressr::progressor(along = focus.n)
        future_sapply(focus.n, 
                      .get.gwGraph,
                      simplify = "array",
                      future.globals = c("kList", "pr"),
                      USE.NAMES = TRUE)
    }
    
    temp <- graphProg_FUN(focus.n = focus.n)
    
    # Set an attribute with all spot names 
    attr(temp, "spot_names") <- attr(kList, "spot_names")
    
    message("step 4/6: DONE!!")
    
    return(temp)
}
