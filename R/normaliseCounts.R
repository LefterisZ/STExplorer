#' Normalise Counts in a Spatial Transcriptomics Experiment
#'
#' This function normalises counts in a spatial transcriptomics experiment
#' using specified methods.
#'
#' @param msfe An object containing spatial transcriptomics experiment data
#' (MetaSpatialFeatureExperiment).
#' @param method The normalisation method to be applied. Currently, only "log"
#' transformation is supported.
#' @param ... Additional arguments to be passed to the normalisation function.
#'
#' @return A normalised SpatialFeatureExperiment object.
#'
#' @details The function normalises counts in a spatial transcriptomics
#' experiment using specified methods. Currently, only the "log" transformation
#' method is supported. Additional arguments can be passed to the underlying
#' normalisation function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords normalisation SpatialFeatureExperiment log transformation
#'
#' @rdname normaliseCounts
#'
#' @examples
#' \dontrun{
#' # Normalise counts using log transformation
#' # msfe <- # MetaSpatialFeatureExperiment object
#' msfe <- normaliseCounts(msfe, method = "log", base = 2)
#' }

normaliseCounts <- function(msfe,
                            method = c("log"),
                            ...) {
  ## Check arguments
  ## When we will establish the metaSFE object remember to let it inherit the
  ## SpatialFeatureExperiment class.
  # stopifnot(is(msfe, "SpatialFeatureExperiment"))

  ## Check valid method argument
  if (missing(method)) {
    method <- "log"
  } else {
    method <- match.arg(method)
  }

  if (method == "log") {
    msfe_int <- lapply(msfe, scater::logNormCounts, ...)
  }

  return(msfe_int)
}
