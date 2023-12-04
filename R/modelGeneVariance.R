#' Model Gene Variance
#'
#' This function models gene variance for each gene in a
#' SpatialFeatureExperiment.
#'
#' @param msfe A SpatialFeatureExperiment object.
#' @param sample_id Either a logical vector indicating the samples to include
#' (if TRUE, all samples are included), or a character vector specifying the
#' sample IDs to include. It is suggested to use a character vector to specify
#' a specific sample.
#' @param method A character string specifying the method for modeling gene
#' variance.
#'   Possible values: "Var" (default), "VarPoisson", "VarSpikes", "CV2",
#'   "CV2Spikes".
#'   \itemize{
#'   \item - "Var": Model gene variance using scran::modelGeneVar.
#'   \item - "VarPoisson": Model gene variance using
#'   scran::modelGeneVarByPoisson.
#'   \item - "VarSpikes": Model gene variance using
#'   scran::modelGeneVarWithSpikes.
#'   \item - "CV2": Model gene coefficient of variation squared using
#'   scran::modelGeneCV2.
#'   \item - "CV2Spikes": Model gene coefficient of variation squared using
#'   scran::modelGeneCV2WithSpikes.
#'   }
#' @param ... Additional arguments to be passed to the chosen method.
#'
#' @return A list of modeled gene variances for each gene.
#'
#' @examples
#' \dontrun{
#' # Model gene variance using the default method ("Var").
#' modelGeneVariance(msfe)
#'
#' # Model gene variance using the "VarPoisson" method.
#' modelGeneVariance(msfe, method = "VarPoisson")
#' }
#'
#' @seealso
#' \code{\link[scran]{modelGeneVar}},
#' \code{\link[scran]{modelGeneVarByPoisson}},
#' \code{\link[scran]{modelGeneVarWithSpikes}},
#' \code{\link[scran]{modelGeneCV2}},
#' \code{\link[scran]{modelGeneCV2WithSpikes}}
#'
#' @references
#' McCarthy DJ, Campbell KR, Lun ATL, Wills QF (2017).
#' Scater: pre-processing, quality control, normalization and visualization of
#' single-cell RNA-seq data in R. Bioinformatics, 33(8), 1179–1186.
#'
#' @importFrom scran modelGeneVar modelGeneVarByPoisson modelGeneVarWithSpikes
#' @importFrom scran modelGeneCV2 modelGeneCV2WithSpikes
#'
#' @export
modelGeneVariance <- function(msfe,
                              sample_id,
                              method = c("Var", "VarPoisson", "VarSpikes",
                                         "CV2", "CV2Spikes"),
                              ...) {
  ## Check arguments
  # stopifnot(is(msfe, "SpatialFeatureExperiment"))
  method <- match.arg(method)

  ## Select samples
  ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)

  ## Model Gene Variance
  var_int <- switch(method,
    Var = lapply(msfe[ids], scran::modelGeneVar, ...),
    VarPoisson = lapply(msfe[ids], scran::modelGeneVarByPoisson, ...),
    VarSpikes = lapply(msfe[ids], scran::modelGeneVarWithSpikes, ...),
    CV2 = lapply(msfe[ids], scran::modelGeneCV2, ...),
    CV2Spikes = lapply(msfe[ids], scran::modelGeneCV2WithSpikes, ...)
  )

  return(var_int)
}

#' Internal Function: Get SpatialFeatureExperiment Sample IDs
#'
#' This internal function retrieves the sample IDs based on the provided
#' criteria.
#'
#' @param msfe A SpatialFeatureExperiment object.
#' @param sample_id Either a logical vector indicating the samples to include
#' (if TRUE, all samples are included), or a character vector specifying the
#' sample IDs to include.
#'
#' @return A character vector of selected sample IDs.
#'
#' @examples
#' \dontrun{
#' # Get all sample IDs in the SpatialFeatureExperiment.
#' .int_getMSFEsmplID(msfe, TRUE)
#'
#' # Get sample IDs based on a logical vector.
#' .int_getMSFEsmplID(msfe, sample_id = c(TRUE, FALSE, TRUE, TRUE))
#'
#' # Get sample IDs based on a character vector.
#' .int_getMSFEsmplID(msfe, sample_id = c("Sample1", "Sample3"))
#' }
#'
.int_getMSFEsmplID <- function(msfe, sample_id) {
  ## Select samples
  if (isTRUE(sample_id)) {
    ids <- names(msfe@sfe_data)
  } else if (is.character(sample_id)) {
    ids <- sample_id
  }

  return(ids)
}
