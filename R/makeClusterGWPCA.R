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
#' @rdname makeClusterGWPCA
#'
#' @examples
#' # Create a parallel cluster for GWPCA computation
#' cluster <- makeClusterGWPCA(spec = 4, type = "PSOCK")
#'
#' @export
makeClusterGWPCA <- function(spec, type = "FORK", ...) {
  ## Get number of clusters
  if (isEmpty(spec) && type == "FORK") {
    spec <- parallelly::availableCores() - 1
  }

  ## Set clusters
  parallel::makeCluster(spec, type = type, ...)

}
