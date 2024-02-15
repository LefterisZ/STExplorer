#' Calculate Per-Gene Logarithmic Means in a Spatial Transcriptomics Experiment
#'
#' This function calculates per-gene logarithmic means in a spatial
#' transcriptomics experiment for each unique sample.
#'
#' @param msfe A MetaSpatialFeatureExperiment object containing spatial
#' transcriptomics experiment data.
#' @param ... Additional arguments. Not used currently
#'
#' @return An updated MetaSpatialFeatureExperiment with a column in rowData
#' for per-gene logarithmic means for each unique sample.
#'
#' @details The function iterates through each unique sample in the
#' MetaSpatialFeatureExperiment and calculates per-gene logarithmic means.
#' Additional arguments can be passed to the underlying calculation function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics experiment, per-gene log mean,
#' MetaSpatialFeatureExperiment
#'
#' @rdname perGeneLogMean
#'
#' @examples
#' \dontrun{
#' # Calculate per-gene logarithmic means in a MetaSpatialFeatureExperiment
#' # msfe <- # A MetaSpatialFeatureExperiment object
#' msfe <- perGeneLogMean(msfe)
#' }
#' @export
perGeneLogMean <- function(msfe,
                           ...) {
  ## Check arguments
  ## When we will establish the metaSFE object remember to let it inherit the
  ## SpatialFeatureExperiment class.
  # stopifnot(is(msfe, "SpatialFeatureExperiment"))

  ## Get the number of unique samples
  ## - Kept here for legacy reasons.
  ## - Will be removed once we transition to metaSFE from the start.
  ids <- names(msfe)

  ## Initialise list
  msfe_int <- list()

  ## Perform calculations
  for (id in ids) {
    msfe_int[[id]] <- .int_perGeneLogMean(id = id, sfe = msfe[[id]])
  }

  return(msfe_int)
}


#' Internal Function: Calculate Per-Gene Logarithmic Mean for a Specific Sample
#' in a MetaSpatialFeatureExperiment object
#'
#' This internal function calculates the per-gene logarithmic mean for a
#' specific sample in a MetaSpatialFeatureExperiment object.
#'
#' @param id The unique identifier of the sample.
#' @param sfe The SpatialFeatureExperiment object for the specified sample.
#'
#' @return The SpatialFeatureExperiment object with per-gene logarithmic means
#' calculated for the specified sample.
#'
#' @details The function calculates the per-gene logarithmic mean for a
#' specific sample in a SpatialFeatureExperiment. It considers the log counts
#' and the number of locations a gene is expressed in. The resulting log mean
#' is saved in the rowData of the SpatialFeatureExperiment object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics experiment, per-gene log mean,
#' internal function, SpatialFeatureExperiment
#'
#' @rdname dot-int_perGeneLogMean
#'
#' @aliases .int_perGeneLogMean
#'
.int_perGeneLogMean <- function(id, sfe) {
  ## Set column names to be used
  colnameLogM <- "s_logMean"
  colnameNLoc <- paste0(id, ".nLocations")
  ## Calculate rowSums only for specific sample
  rowSum <- rowSums(assay(sfe, "logcounts"))
  ## Get number of locations
  nLocations <- rowData(sfe)[[colnameNLoc]]
  ## Calculate the mean of log counts over the locations a gene is expressed
  sLogMean <- rowSum / nLocations
  ## We must account for the non-expressed genes which will have a mean of NaN
  sLogMean <- sLogMean %>%
    data.frame(sLogMean = .) %>%
    dplyr::mutate(sLogMean = dplyr::if_else(is.na(sLogMean), 0, sLogMean))

  ## Save the log means
  rowData(sfe)[[colnameLogM]] <- sLogMean
  ## Return
  return(sfe)
}
