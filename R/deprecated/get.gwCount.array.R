#' @name get.gwCount.array
#' 
#' @description A function to produce a 3D array with dimensions s x g x s. 
#'              Where s = spots (locations) and g = genes. Takes as input an 
#'              s x g matrix of gene counts, an s x s distance matrix with 
#'              weighted distances between locations. Returns as output a 3D
#'              array with gene expression counts weighted for each location (or
#'              for a number of selected locations).
#' 
#' @param obs A matrix containing the observation data, usually gene expression 
#'            data where columns are the genes and rows are the locations.
#' @param wdmat An n x n matrix (n = number of locations) of weighted distances 
#'              between all locations.
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
#' @return A 3D array with dims = [s, g, focus.n]
#' 
#' @export


get.gwCount.array <- function(obs, wdmat, focus.n, strategy = "multicore", 
                              workers = 5){
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:nrow(obs) # indexes of locations to use (z-axis)
        names(focus.n) <- rownames(obs) # give spot names to the indexes
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of location indexes was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
        names(focus.n) <- dimnames(obs)[[1]][focus.n] # give names
    } else if(is.vector(focus.n) & is.character(focus.n)) {
        message("A selection of location names was provided...")
        message("Locations with names: ", paste(focus.n, collpse = " "))
        focus.n <- match(focus.n, dimnames(obs)[[1]]) # find the indexes
        names(focus.n) <- dimnames(obs)[[1]][focus.n] # give names
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values OR character values only. 
             Please make sure you provide location indexes only.")
    }
    
    # # Weight expression data --> sequential and slow
    # temp <- sapply(focus.n, .get.gw.counts,
    #                obs = obs, wdmat = wdmat,
    #                simplify = "array")
    
    # Weight expression data in parallel with future package
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
    
    wCountProg_FUN <- function(focus.n, wdmat){
        pr <- progressr::progressor(along = focus.n)
        future_sapply(focus.n, .get.gw.counts,
                      wdmat = wdmat,
                      simplify = "array",
                      future.globals = c("obs", "pr", "temp.list"),
                      future.chunk.size = dim(wdmat)[1]/workers,
                      USE.NAMES = TRUE)
    }

    temp <- wCountProg_FUN(focus.n = focus.n,
                           wdmat = wdmat)

    plan(strategy = "sequential")
    
    ## add dimnames --> spot names as rows and columns
    # dimnames(temp)[[1]] <- rownames(obs)
    # dimnames(temp)[[2]] <- colnames(obs)

    message("step 1/6: DONE!!")
    
    return(temp)
    
}
