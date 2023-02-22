#' @name get.gwGraph.list
#' 
#' @description a function to wrap get.gwGraph.table function for parallel 
#'              processing using foreach, progressr and doFuture packages.
#' @param .focus.n a character vector of indexes for the locations in focus.
#' @param .obs the gene expression matrix
#' @param .wdmat the geographically weighted physical distance matrix
#' @param .method the method to be used for the statistical distance. Default is 
#'               "euclidean".
#' @param .k k for the k-NN. Keep in mind that a self neighbour is calculated and
#'          then dropped. So, if you want 7 neighbours you must set k = 8.
#' @param .combine combine reciprocal edges. This means that A --> B is 
#'                 considered the same as B --> A. Default is FALSE.
#' @param .strategy the future plan strategy. More info at future::plan
#' @param .workers the number of cores to be used. More info at future::plan
#' @param verbose default TRUE. Show progress bar.
#' @param ... passes arguments to the Rfast::Dist function
#'

get.gwGraph.list <- function(.focus.n,
                             .obs,
                             .wdmat,
                             .method,
                             .k,
                             .verbose,
                             .combine,
                             ...) {
  
  ## Set progress to be visible
  if (.verbose) {
    ## Run verbose
    progrFUN(..focus.n = .focus.n,
             ..obs = .obs,
             ..wdmat = .wdmat,
             ..method = .method,
             ..k = .k,
             ..combine = .combine,
             ...)
    
  } else {
    foreach::foreach(i = .focus.n) %dopar% {  
      ## Create the graph table
      get.gwGraph.table(..focus = i, 
                        ..obs = .obs,
                        ..wdmat = .wdmat,
                        ..method = .method,
                        ..k = .k,
                        ..combine = .combine,
                        ...)
    }
  }
}

# ------------------------------------------- #
# ------ internally defined functions ------- #

progrFUN <- function(..focus.n, 
                     ..obs,
                     ..wdmat,
                     ..method,
                     ..k,
                     ..combine,
                     ...) {
  ## Set progress handler
  pr <- progressor(along = ..focus.n)
  ## Run the foreach
  foreach::foreach(i = ..focus.n) %dopar% { 
    ## Fetch progress
    pr()
    ## Create the graph table
    get.gwGraph.table(..focus = i, 
                      ..obs,
                      ..wdmat,
                      ..method,
                      ..k,
                      ..combine,
                      ...)
  }
}

