#'  @name get.gwGraphs
#'  
#'  @description A function that will generate a list of GW graphs.
#'  
#' @param focus.n a character vector of indexes for the locations in focus. 
#'                Default is NULL and all locations are considered.
#' @param obs the gene expression matrix
#' @param wdmat the geographically weighted physical distance matrix
#' @param method the method to be used for the statistical distance. Default is 
#'               "euclidean".
#' @param k k for the k-NN. Keep in mind that a self neighbour is calculated and
#'          then dropped. So, if you want 7 neighbours you must set k = 8.
#' @param combine combine reciprocal edges. This means that A --> B is 
#'                 considered the same as B --> A. Default is FALSE.
#' @param strategy the future plan strategy. More info at future::plan
#' @param workers the number of cores to be used. More info at future::plan
#' @param verbose default TRUE. Show progress bar.
#' @param ... passes arguments to the Rfast::Dist function
#' 
#' @export

get.gwGraphs <- function(focus.n = NULL,
                         obs,
                         wdmat,
                         method = "euclidean",
                         k = 7,
                         combine = FALSE,
                         strategy = "sequential",
                         workers = 1,
                         verbose = TRUE,
                         ...){
  
  # Set handlers for a verbose approach ---
  if (verbose) {
    set.verbose(.global = TRUE,
                .handlers = c("progress", "beepr"))
  }
  
  # Sort the observations table and generate the dictionaries ----
  dictionary <- get.dictionary(obs)
  
  # Change the row and column names based on the dictionary ----
  obs <- add.indexes(obs = obs, dict_list = dictionary, row = TRUE, col = TRUE)
  
  # Sort the weighted distances matrix by spot name ----
  ## this is sorted at the wdmat generation code snippet. We sort the polygons 
  ## object by the "Barcode" column and then calculate physical distances and 
  ## apply gw.weights. The resulting wdmat matrix is without row and column names. 
  ## So, we add row and column names ranging from 1:nrow(wdmat) and 1:ncol(wdmat) 
  ## like we did for the polygons object's rownames because we know that the way 
  ## gw.dist and gw.weight functions work, they keep the given order of locations.
  
  # Check indexes of locations ----
  focus.n <- check.focus.n(.focus.n = focus.n, .dict = dictionary)
  
  # Set parallelisation if an appropriate strategy is supplied  ----
  ## Fetch the user's current backend ----
  oplan <- plan()
  
  ## Set parallel or sequential backend ----
  check.parallel(.strategy = strategy, .workers = workers)
  
  ## doFuture is using a parallel-safe RNG method as per this:
  ## https://github.com/tidymodels/tune/issues/377 therefore we silence the 
  ## warning.
  if (getDoParWorkers() > 1) {
    rlang::local_options(doFuture.rng.onMisuse = "ignore")
  }
  
  # Prepare the user for long waiting times ----
  message("GW graphs to be generated: ", length(focus.n))
  message("Please be patient. This will take a while...")
  mins <- round((length(focus.n) * 3.2) / 60, 3)
  message("E.T.A. is approximatelly: ", mins, " minutes")
  message("Sit back and relax!! :)")
  ## fetch start time
  s <- Sys.time()
  
  # Make the graphs list ----
  graphs_list <- get.gwGraph.list(.focus.n = focus.n,
                                  .obs = obs,
                                  .wdmat = wdmat,
                                  .method = method,
                                  .k = k,
                                  .verbose = verbose,
                                  .combine = combine)
  
  ## set the 'focus' attribute as names for the list's objects 
  graphs_list <- set.listNames(list = graphs_list)
  
  ## fetch end time
  e <- Sys.time()
  ## return diference
  message("Generating GW graphs is done!")
  message("Time elapsed: ", round(difftime(e, s, units = "mins"), 3), " mins")
  
  
  # Reset the user's backend ----
  ## A future plan should never be called inside a function. But, because here we
  ## need to do so, we first save the user's backend and then we set it back.
  
  ## Stop cluster ----
  if (strategy == "cluster" | inherits(strategy, "ClusterFuture")) {
    parallel::stopCluster(workers)
  }
  ## Reset plan ----
  plan(oplan)
  plan()
  
  ## Return the list ----
  return(graphs_list)
}
