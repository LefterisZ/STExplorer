#' @name get.gwDist.array
#' 
#' @description A wrapper for base R dist function. Gives the ability to return
#'              a 3D array with dimensions s x s x s. Where s = spots (locations). 
#'              Takes as input an s x g (genes) matrix of spatially weighted 
#'              gene counts. Returns as output a 3D array with statistical 
#'              distances between locations based on the weighted gene expression
#'              calculated by the get.gwCount.array() function.
#'              
#' @param gwCounts The input data, must be a 3D array of gw-counts as generated
#'                 by the get.gwCount.array() function.
#' @param focus.n Numeric or character. The indexes (numeric) or the names 
#'                (character) of the locations you want in focus. It can be a 
#'                vector of indexes from locations that are of interest. Default
#'                behaviour is for focus.n to be missing. This will result to  
#'                all locations being considered.
#' @param method The distance measure to be used. This must be one of "euclidean", 
#'               "maximum", "manhattan", "canberra", "binary" or "minkowski". 
#'               Any unambiguous substring can be given. For more info, look at 
#'               the base R dist function. Default is euclidean.
#' @param p The power of the Minkowski distance. For more info, look at the base 
#'          R dist function. Default is 2 which is the euclidean distance.
#' 
#' @param strategy a future plan strategy for parallelisation. The defaults is 
#'                 multicore because it doesn't need to load the environment on 
#'                 each worker and thus is more suitable for large datasets. For
#'                 more info look at future package.
#' @param workers the number of cores to use for parallelisation.
#' 
#' @return a 3D array with dims = [s, s, focus.n]
#' 
#'
get.gwDist.array <- function(gwCounts, focus.n, method = "euclidean", p = 2, 
                             strategy = "multicore", workers = 5){
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:dim(gwCounts)[3] # indexes of locations to use (z-axis)
        names(focus.n) <- rownames(gwCounts) # give spot names to the indexes
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of location indexes was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
        names(focus.n) <- dimnames(gwCounts)[[1]][focus.n] # give names
    } else if(is.vector(focus.n) & is.character(focus.n)) {
        message("A selection of location names was provided...")
        message("Locations with names: ", paste(focus.n, collpse = " "))
        focus.n <- match(focus.n, dimnames(gwCounts)[[3]]) # find the indexes
        names(focus.n) <- dimnames(gwCounts)[[1]][focus.n] # give names
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values OR character values only. 
             Please make sure you provide location indexes only.")
    }
    
    # Find distances without parallelisation
    # temp <- sapply(focus.n,
    #                function(X, method, p){
    #                    message("Calculating location: ", X)
    #                    dist(gwCounts[,,X], method, p) %>% as.matrix()
    #                },
    #                method = method,
    #                p = p,
    #                simplify = "array")
    
    
    # Find the distances with future parallelisation
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
    
    # obj.s <- object.size(gwCounts)*2
    # options("future.globals.maxSize" = obj.s)
    
    message("Calculating the distance matrices...")
    message("This might take some time depending on the number of spots and 
            genes you included.")

    distProg_FUN <- function(focus.n, method, p){
        pr <- progressr::progressor(along = focus.n)
        future_sapply(focus.n,
                      .dist.wrap,
                      method = method,
                      p = p,
                      # gwCounts = gwCounts,
                      simplify = "array",
                      future.globals = c("gwCounts", "pr"),
                      USE.NAMES = TRUE)
    }


    temp <- distProg_FUN(focus.n = focus.n,
                         method = method,
                         p = p)
    
    plan(strategy = "sequential")
    
    ## add dimnames --> spot names as rows and columns
    dimnames(temp)[[1]] <- rownames(gwCounts)
    dimnames(temp)[[2]] <- rownames(gwCounts)
   
    message("step 2/6: DONE!!") 
    
    return(temp)
}

