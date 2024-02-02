#' Filter In-Tissue Spots in Spatial Transcriptomics Data
#'
#' This function filters in-tissue spots in spatial transcriptomics data,
#' providing users with the flexibility to specify sample/image identifiers.
#' The sample_id parameter accepts TRUE (equivalent to all samples/images) or
#' NULL (specifying the first available entry).
#'
#' The function internally utilises the .int_filterInTissue function to filter
#' in-tissue spots for each specified sample. The results are then merged into a
#' single SpatialFeatureExperiment object using the Reduce function.
#'
#' Additionally, the function ensures data integrity by cleaning up duplicate
#' entries in the metadata slot with the .int_cleanMetaData function.
#'
#' This filtering step is crucial for focusing analyses on relevant in-tissue
#' spots in spatial transcriptomics experiments.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param sample_id A character string, TRUE, or NULL specifying sample/image
#'        identifier(s).
#'
#' @return A \code{SpatialFeatureExperiment} object with in-tissue spots
#' filtered.
#'
#' @seealso \code{\link{.int_filterInTissue}}, \code{\link{.int_cleanMetaData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics, in-tissue spots, data filtering
#'
#' @rdname filterInTissue
#'
#' @examples
#' \dontrun{
#' sfe <- filterInTissue(sfe)
#' }
#'
#'
#' @export
filterInTissue <- function(sfe, sample_id = TRUE) {
  ## Get a vector of sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Create a list of filtered SpatialFeatureExperiment objects
  sfe_list <- lapply(ids, .int_filterInTissue, sfe = sfe)

  ## Merge the list of SpatialFeatureExperiment objects into one
  sfeOut <- Reduce(function(x, y) cbind(x, y), sfe_list)

  ## Clean up duplicate entries in the metadata slot
  sfeOut <- .int_cleanMetaData(sfeOut = sfeOut)

  return(sfeOut)

}


#' Internal Function: Filter In-Tissue Spots for a Specific Sample
#'
#' This internal function filters in-tissue spots for a specific sample in
#' spatial transcriptomics data. It takes a character string (\code{id})
#' specifying the sample/image identifier and a \code{SpatialFeatureExperiment}
#' object (\code{sfe}). The logic is implemented using the .int_filterInTissue
#' function, and the resulting \code{SpatialFeatureExperiment} object contains
#' only in-tissue spots for the specified sample.
#'
#' @param id A character string specifying the sample/image identifier.
#' @param sfe A \code{SpatialFeatureExperiment} object.
#'
#' @return A \code{SpatialFeatureExperiment} object with in-tissue spots
#' filtered for the specified sample.
#'
#' @details
#' This internal function efficiently filters in-tissue spots for a specific
#' sample in spatial transcriptomics data. It takes a character string
#' (\code{id}) specifying the sample/image identifier and a
#' \code{SpatialFeatureExperiment} object (\code{sfe}). The logic is based on
#' the .int_filterInTissue function, and the resulting
#' \code{SpatialFeatureExperiment} object contains only in-tissue spots for the
#' specified sample.
#'
#'
#' @keywords internal function, spatial transcriptomics, in-tissue spots
#'
.int_filterInTissue <- function(id, sfe) {
  sfe_int <- sfe[, colData(sfe)$sample_id == id]
  in_tissue <- Matrix::which(colData(sfe_int)$in_tissue == TRUE)
  sfe_int <- sfe_int[, in_tissue]
  return(sfe_int)
}


#' Internal Function: Clean Metadata Entries
#'
#' This internal function cleans metadata entries in a SpatialFeatureExperiment
#' object. It ensures that duplicate entries in the metadata slot are removed,
#' preserving the original names. The function takes a SpatialFeatureExperiment
#' object (sfeOut) as input and returns the cleaned object.
#'
#' @param sfeOut A \code{SpatialFeatureExperiment} object.
#'
#' @return A \code{SpatialFeatureExperiment} object with cleaned metadata
#' entries.
#'
#' @details
#' This internal function efficiently cleans metadata entries in a
#' SpatialFeatureExperiment object. It takes a SpatialFeatureExperiment object
#' (sfeOut) as input and ensures that duplicate entries in the metadata slot
#' are removed. The function maintains the original names while making them
#' unique.
#'
#' @keywords internal function, spatial transcriptomics, metadata
#'
.int_cleanMetaData <- function(sfeOut) {
  ## Clean up duplicate entries in the metadata slot
  ## Get names
  mdt_names <- unique(names(metadata(sfeOut)))
  ## Make names unique
  names(metadata(sfeOut)) <- make.unique(names(metadata(sfeOut)))
  ## Keep slots with the original names
  metadata(sfeOut) <- metadata(sfeOut)[mdt_names]

  return(sfeOut)
}
