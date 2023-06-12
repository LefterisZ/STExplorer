#' @name set.parallel
#' 
#' @description 
#' this function is used internally to check and set a backend for gwc and gwpca
#' 
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
    if (getDoParRegistered()) {
      message(sprintf('Currently using %s [%s]\n', dpN, dpV))
    }
    
    ### check if registered
    message(sprintf('%s backend is registered\n',
                    if (getDoParRegistered()) 'A' else 'No'))
    
    ### check the number of execution workers
    message(sprintf('Running with %d worker(s)\n', getDoParWorkers()))
    
  } else {
    ## Run foreach sequentially
    message("Running sequentially.\n")
    foreach::registerDoSEQ()
  }
}
