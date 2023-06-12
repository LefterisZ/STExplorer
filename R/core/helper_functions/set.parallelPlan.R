#' @name set.parallelPlan
#' 
#' @description a function to check the parallelisation parameters and set up 
#'              the environment with the future package.
#' 
#' 

set.parallelPlan <- function(..strategy, ..workers, ...){
  
  plans <- c("multicore", "multisession", "cluster", "sequential", "clustAuto")
  
  ## Find number of available cores
  availCores <- availableCores()
  
  ## Check if strategy is list
  if (inherits(..strategy, c("Future", "environment"))) {
    strat <- ..strategy[["calls"]][[1]][[1]][[3]]
  } else if (..strategy %in% plans) {
    strat <- ..strategy
  } else {
    stop("Please provide a valid future plan (see future::plan).")
  }
  
  ## Check if workers is a cluster object
  if (inherits(..workers, c("SOCKcluster", "cluster"))) {
    work <- TRUE
    demandCores <- length(..workers)
  } else {
    work <- FALSE
    demandCores <- ..workers
  }
  
  ## Check if there are enough cores
  if (availCores < 3 & strat != "sequential") {
    message("Number of cores available: ", availCores)
    stop("You don't have enough cores to parallelise. Please select the\n 
    sequential strategy")
    
  } else if (availCores >= 3 & demandCores >= availCores) {
    message("You have less cores available than the number of workers you\n
    asked for. Or you are asking to use all of your cores.")
    message("Number of cores available: ", availCores)
    message("Number of workers asked: ", demandCores)
    stop("Please change the number of workers.")
  }
  
  ## Check for strategy: cluster
  if (strat == "cluster" & work) {
    ## Set up strategy
    plan(strategy = ..strategy, workers = ..workers, ...)
    message("------------------------------------------------------")
    message("Future parallel strategy: ", strat)
    message("Number of cores to be used: ", demandCores, "/", availCores)
    message("------------------------------------------------------")
    
  } else if (strat == "cluster" & !work) {
    message("A cluster object should be provided with future strategy cluster.")
    message("Please provide a valid cluster object for the 'workers' argument")
    stop("The cluster object can be created like this:
            
        my.cl <- parallel::makeCluster(availableCores() - 1, type = 'FORK')\n")
    
    ## If strategy != cluster check that 'workers' is numeric:
  } else if (strat != "cluster" & is.numeric(..workers)) {
    ## Set up strategy
    plan(strategy = ..strategy, workers = ..workers)
    message("------------------------------------------------------")
    message("Future parallel strategy: ", strat)
    message("Number of cores to be used: ", ..workers, "/", availCores)
    message("------------------------------------------------------")
    
  } else {
    stop("Please provide a valid numeric value for 'workers' argument.")
  }
  
  ## Check for auto-assignment of clusters
  if (strat == "clustAuto") {
    ..workers <- parallel::makeCluster(availableCores() - 2, type = 'FORK')
    plan(strategy = "cluster", workers = ..workers, ...)
    demandCores <- availableCores() - 2
    message("Clusters generated automatically with:
              parallel::makeCluster(availableCores() - 2, type = 'FORK')\n")
    message("------------------------------------------------------")
    message("Future parallel strategy: ", strat)
    message("Number of cores to be used: ", demandCores, "/", availCores)
    message("------------------------------------------------------")
  }
  
}
