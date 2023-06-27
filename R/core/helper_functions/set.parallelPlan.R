#' Set parallel plan
#'
#' @name set.parallelPlan
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' It is used internally to check the parallelisation parameters and set up the
#' environment with the future package.
#'
#' @param ..strategy The strategy for parallel processing. It can be one of the
#' following: "multicore", "multisession", "cluster", "sequential",
#' or "clustAuto".
#'
#' @param ..workers The number of workers to be used for parallel processing. If
#' the strategy is set to "cluster", this should be a cluster object. Otherwise,
#' it should be a numeric value specifying the number of workers. For the
#' "clustAuto" strategy, this parameter is ignored.
#'
#' @param ... Further arguments passed to \code{\link{plan}} function.
#'
#' @importFrom future availableCores plan
#'
#' @details
#' This function checks the specified parallelisation parameters and sets up the
#' environment with the future package. It verifies the availability of cores,
#' checks the strategy and workers, and sets the future plan accordingly. If the
#' strategy is set to "cluster", a valid cluster object should be provided as
#' the workers argument. For other strategies, a numeric value specifying the
#' number of workers should be provided. The function provides informative
#' messages about the chosen strategy and the number of cores to be used.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
set.parallelPlan <- function(..strategy, ..workers, ...){

  plans <- c("multicore", "multisession", "cluster", "sequential", "clustAuto")

  ## Find number of available cores
  availCores <- future::availableCores()

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
    future::plan(strategy = ..strategy, workers = ..workers, ...)
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
    future::plan(strategy = ..strategy, workers = ..workers)
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
    future::plan(strategy = "cluster", workers = ..workers, ...)
    demandCores <- availableCores() - 2
    message("Clusters generated automatically with:
              parallel::makeCluster(availableCores() - 2, type = 'FORK')\n")
    message("------------------------------------------------------")
    message("Future parallel strategy: ", strat)
    message("Number of cores to be used: ", demandCores, "/", availCores)
    message("------------------------------------------------------")
  }

}
