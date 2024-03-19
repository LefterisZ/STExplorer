# ---------------------------------------------------------------------------- #
#  ############ FUNCTION ASSOCIATED WITH GWPCA PARALLELISATION ##############
# ---------------------------------------------------------------------------- #
#' Create a Parallel Cluster for GWPCA Computation
#'
#' This function creates a parallel cluster for running the GWPCA
#' (Geographically Weighted Principal Component Analysis) computation in
#' parallel on multiple cores.
#'
#' @param spec The number of clusters to be created. If not provided and type
#' is "FORK", it defaults to availableCores() - 1.
#' @param type The type of cluster to be created. Default is "FORK". Other
#' options may include "PSOCK" or others supported by parallel::makeCluster.
#' @param ... Additional arguments to be passed to parallel::makeCluster.
#'
#' @return A parallel cluster for GWPCA computation.
#'
#' @details The function sets up a parallel cluster for running the GWPCA
#' computation. The number of clusters can be specified, and the default is to
#' use availableCores() - 1 if none is provided.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords parallel computing, cluster, GWPCA
#'
#' @importFrom parallelly availableCores
#'
#' @rdname makeClusterGWPCA
#'
#' @examples
#' # Create a parallel cluster for GWPCA computation
#' cluster <- makeClusterGWPCA(spec = 4, type = "PSOCK")
#'
#' @export
makeClusterGWPCA <- function(spec = NULL, type = "FORK", ...) {
  ## Get number of clusters
  if (is.null(spec) && type == "FORK") {
    spec <- parallelly::availableCores() - 1
  }

  ## Set clusters
  parallel::makeCluster(spec, type = type, ...)

}


#' Set up verbose handlers
#'
#' @name set.verbose
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' It is used internally to set the requirements for verbosity using the
#' progressr package.
#'
#' @keywords internal
#'
#' @param .global Logical value indicating whether to set the global handler for
#' verbosity. If TRUE, the global handler will be set; if FALSE, it will not
#' be set.
#' @param .handlers A character vector specifying the types of handlers to be
#' set for verbosity. Possible values include "progress" and "beepr".
#'
#' @importFrom progressr handlers
#'
#' @details
#' This function sets up the requirements for verbosity through the progressr
#' package. By default, it sets the global handler for verbosity, which affects
#' all verbosity calls. Additionally, it allows specifying specific types of
#' handlers to be used for verbosity, such as "progress" and "beepr". The
#' function helps control the verbosity settings in the code.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
set.verbose <- function(.global = TRUE,
                        .handlers = c("progress", "beepr")){
  ## Set global handler?
  progressr::handlers(global = .global)
  ## Set progress types
  progressr::handlers(.handlers)
}


# ---------------------------------------------------------------------------- #
#  ############### FUNCTIONS ASSOCIATED WITH PARALLELISATION ################
# ---------------------------------------------------------------------------- #
#' Chekc and set backend
#'
#' @name set.parallel
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' It is used internally to check and set a backend for the "gwcSTE" and
#' "gwpcaSTE" functions.
#'
#' @keywords internal
#'
#' @param .strategy The strategy to be used for parallel processing. If not
#' specified or set to "sequential", the function will run sequentially.
#'
#' @param .workers The number of workers to be used for parallel processing.
#' Ignored if the strategy is set to "sequential".
#'
#' @details
#' This function sets the parallel processing backend based on the specified
#' strategy. If the strategy is set to "sequential" or if no strategy is
#' specified, the function will run sequentially using the "registerDoSEQ"
#' backend. If a different strategy is specified, the function will set the
#' future plan and register the appropriate backend using the "registerDoFuture"
#' function from the "doFuture" package. The number of workers can also be
#' specified, which determines the level of parallelism.
#'
#' @importFrom doFuture registerDoFuture
#' @importFrom foreach getDoParVersion getDoParName getDoParRegistered
#' @importFrom foreach getDoParWorkers registerDoSEQ
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
set.parallel <- function(.strategy, .workers) {

  if (!.strategy == "sequential" | !inherits(.strategy, "SequentialFuture")) {
    ## set the future plan
    set.parallelPlan(..strategy = .strategy, ..workers = .workers)

    ## Run foreach with doFuture's %dopar%
    ### register doFuture to use %dopar%
    doFuture::registerDoFuture()

    ### check the name of the currently registered doPar backend
    dpV <- foreach::getDoParVersion()
    dpN <- foreach::getDoParName()
    if (foreach::getDoParRegistered()) {
      message(sprintf('Currently using %s [%s]\n', dpN, dpV))
    }

    ### check if registered
    message(sprintf('%s backend is registered\n',
                    if (foreach::getDoParRegistered()) 'A' else 'No'))

    ### check the number of execution workers
    message(sprintf('Running with %d worker(s)\n', foreach::getDoParWorkers()))

  } else {
    ## Run foreach sequentially
    message("Running sequentially.\n")
    foreach::registerDoSEQ()
  }
}


#' Set parallel plan
#'
#' @name set.parallelPlan
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' It is used internally to check the parallelisation parameters and set up the
#' environment with the future package.
#'
#' @keywords internal
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
    stop("You don't have enough cores to parallelise. Please select the ",
         "sequential strategy")

  } else if (availCores >= 3 & demandCores >= availCores) {
    message("You have less cores available than the number of workers you ",
            "asked for. Or you are asking to use all of your cores.")
    message("Number of cores available: ", availCores)
    message("Number of workers asked: ", demandCores)
    stop("Please change the number of workers.")
  }

  ## Check for strategy: sequential
  if (strat == "sequential") {
    ## Set up strategy
    future::plan(strategy = ..strategy, ...)
  }

  ## Check for strategy: cluster
  if (strat == "cluster" & work) {
    ## Set up strategy
    future::plan(strategy = ..strategy, workers = ..workers, ...)

  } else if (strat == "cluster" & !work) {
    message("A cluster object should be provided with future strategy cluster.")
    message("Please provide a valid cluster object for the 'workers' argument")
    stop("The cluster object can be created like this:\n",
         "\n",
         "        my.cl <- parallel::makeCluster(",
         "parallelly::availableCores() - 1, type = 'FORK')\n")

    ## If strategy != cluster check that 'workers' is numeric:
  } else if (strat != "cluster" & is.numeric(..workers)) {
    ## Set up strategy
    future::plan(strategy = ..strategy, workers = ..workers)

  } else {
    stop("Please provide a valid numeric value for 'workers' argument.")
  }

  ## Check for auto-assignment of clusters
  if (strat == "clustAuto") {
    ..workers <- parallel::makeCluster(availableCores() - 2, type = 'FORK')
    future::plan(strategy = "cluster", workers = ..workers, ...)
    demandCores <- availableCores() - 2
    message("Clusters generated automatically with:\n",
            "parallel::makeCluster(",
            "parallelly::availableCores() - 1, type = 'FORK')\n")
  }

  ## Print out set-up message
  message("\n------------------------------------------------------")
  message("Future parallel strategy: ", strat)
  message("Number of cores to be used: ", demandCores, "/", availCores)
  message("------------------------------------------------------")

}
