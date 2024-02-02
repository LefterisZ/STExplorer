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

# The internal function `.int_getMSFEsmplID` can be found in the
# `class_metaSFE.R` file.


#' Get Top High Variable Genes from a Spatial Transcriptomics Experiment
#'
#' This function retrieves the top high variance genes from a spatial
#' transcriptomics experiment based on a provided DataFrame of variance
#' modelling statistics.
#'
#' @param stats A \code{DataFrame} of variance modelling statistics with one
#' row per gene. Alternatively, a SpatialFeatureExperiment object, in which
#' case it is supplied to modelGeneVar to generate the required DataFrame.
#' @param sample_id TRUE selects all samples. Alternatively, provide a
#' character vector with the sample IDs for filtering.
#' @param ... Additional arguments to be passed to
#' \code{\link[scran]{getTopHVGs}}.
#'
#' @return A list of top high variance genes for each sample, as determined by
#' the provided DE analysis results.
#'
#' @details The function uses a DE analysis result list (`stats`) to extract
#' the top high variance genes for each sample. It allows for optional
#' filtering based on sample IDs. Additional arguments can be passed to
#' the underlying \code{\link[scran]{getTopHVGs}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords high variance genes, spatial transcriptomics experiment
#'
#' @seealso [getTopHVGs()]
#'
#' @rdname getTopHighVarGenes
#'
#' @examples
#' \dontrun{
#' # Get top high variance genes from a spatial transcriptomics experiment
#' # 'stats' a DataFrame of variance modelling statistics
#' top_hvg <- getTopHighVarGenes(stats, sample_id = TRUE)
#' }
#'
#' @export
getTopHighVarGenes <- function(stats, sample_id = TRUE, ...) {
  ## Check arguments
  # stopifnot(is(msfe, "SpatialFeatureExperiment"))
  stopifnot(is.list(stats))

  ## Select samples
  ids <- .int_getMSFEsmplID(msfe = stats, sample_id = sample_id)

  ## Get top high variable genes
  hvgs_int <- lapply(stats[ids], scran::getTopHVGs, ...)

  return(hvgs_int)
}

