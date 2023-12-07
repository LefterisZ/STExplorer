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
  ids <- .int_getMSFEsmplID(list = stats, sample_id = sample_id)

  ## Get top high variable genes
  hvgs_int <- lapply(stats[ids], scran::getTopHVGs, ...)

  return(hvgs_int)
}
