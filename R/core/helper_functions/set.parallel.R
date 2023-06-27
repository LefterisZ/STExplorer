#' Chekc and set backend
#'
#' @name set.parallel
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' It is used internally to check and set a backend for the "gwcSTE" and
#' "gwpcaSTE" functions.
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
