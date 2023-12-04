#' Define the MetaSpatialFeatureExperiment class using S4 classes and methods
#'
#' This class represents a collection of SpatialFeatureExperiment (SFE) objects
#' for spatial transcriptomics data. Each SFE object is stored in a slot within
#' the MetaSpatialFeatureExperiment object, and the slot names correspond to
#' the experiment/sample names.
#'
#' @export
setClass(
  "MetaSpatialFeatureExperiment",
  representation(
    sfe_data = "list",
    sample_ids = "list"
  )
)

#' Constructor function for MetaSpatialFeatureExperiment class
#'
#' @param sfe_data A list of SpatialFeatureExperiment objects.
#' @param sample_ids A list of sample IDs matching the sample IDs of the stored
#' SpatialFeatureExperiment (SFE) objects.
#'
#' @return An instance of the MetaSpatialFeatureExperiment class.
#'
#' @examples
#' ## Create a MetaSpatialFeatureExperiment object
#' msfe <- MetaSpatialFeatureExperiment()
#'
#' ## Add an SFE object to the MetaSpatialFeatureExperiment
#' ## create or load your SpatialFeatureExperiment object
#'
#' # sfe <- ...
#'
#' ## a SpatialFeatureExperiment object encompassing a single sample named, in
#' ##  this instance, 'sample_1'.
#'
#' # msfe <- addSFE(msfe, sfe)
#'
#' ## Access an SFE object by name
#' # sfe_sample_1 <- getSFE(msfe, "sample_1")
#'
#' @export
MetaSpatialFeatureExperiment <- function(sfe_data = list(),
                                         sample_ids = list()) {
  new("MetaSpatialFeatureExperiment",
      sfe_data = sfe_data,
      sample_ids = sample_ids)
}

#' Method to add a SpatialFeatureExperiment to the MetaSpatialFeatureExperiment
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param sfe An instance of the SpatialFeatureExperiment class to add.
#' @param sample_id Currently only accepts TRUE. Internally performs automatic
#' selection of the sample ID in the SFE object. Assumes that the SFE object
#' includes only one sample inside. If multiple samples are stored in the same
#' SFE object, please use the \code{addMultipleSFE} instead.
#'
#' @return An updated MetaSpatialFeatureExperiment object with the added SFE.
#'
#' @export
addSFE <- function(x, sfe, sample_id = TRUE) {
  id <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)
  x@sfe_data[[id]] <- sfe
  x
}

#' Method to access a SpatialFeatureExperiment by name
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param name Name of the experiment/sample for the SFE to access.
#'
#' @return The SFE object associated with the provided name.
#'
#' @export
getSFE <- function(x, name) {
  x@sfe_data[[name]]
}

#' Method to add an SFE object with multiple samples to the
#' MetaSpatialFeatureExperiment
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param sfe An instance of the SpatialFeatureExperiment class with multiple
#' samples.
#'
#' @return An updated MetaSpatialFeatureExperiment object with the added SFEs.
#'
#' @export
addMultipleSFE <- function(x, sfe) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = TRUE)

  for (id in ids) {
    ## Subset per sample
    sfe_int <- sfe[, colData(sfe)$sample_id == id]

    ## Add each subset to the MetaSpatialFeatureExperiment object
    x <- addSFE(x, sfe_int, sample_id = id)
  }

  x
}

#' Method to retrieve multiple SFE objects as a single SFE object
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param sample_ids The sample IDs from the SFE objects to be combined and
#' retrieved as one SFE object
#'
#' @return An SFE object with multiple samples inside.
#'
#' @export
getMultipleSFE <- function(x, sample_ids) {
  ## Check provided sample IDs are correct

  ## Combine into one SFE object


  return(sfe)
}

#' Method to retrieve sample IDs from MetaSpatialFeatureExperiment
#'
#' This method retrieves the sample IDs from a MetaSpatialFeatureExperiment object.
#' The sample IDs are stored in the @sample_ids slot of the MetaSpatialFeatureExperiment.
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#'
#' @return A character vector containing the sample IDs.
#'
#' @export
#' @examples
#' # Create a MetaSpatialFeatureExperiment object
#' msfe <- MetaSpatialFeatureExperiment()
#'
#' # Create or load a list of SpatialFeatureExperiment objects
#' # Replace sfe1, sfe2, sfe3 with your actual SFE objects
#' # sfe_list <- list(sfe1, sfe2, sfe3)
#'
#' # Add multiple SFE objects to the MetaSpatialFeatureExperiment
#' # msfe <- addMultipleSFEs(msfe, sfe_list)
#'
#' # Retrieve sample IDs
#' sample_ids <- getSampleIDs(msfe)
#'
#' @seealso \code{\link{MetaSpatialFeatureExperiment}}, \code{\link{addMultipleSFEs}}
#' @keywords internal
getSampleIDs <- function(x) {
  unlist(x@sample_ids)
}

