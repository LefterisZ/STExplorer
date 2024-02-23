#' Apply Quality Control Thresholds to SpatialFeatureExperiment
#'
#' This function applies quality control (QC) thresholds to a
#' \code{SpatialFeatureExperiment} object. It filters out features (genes or
#' other genomic elements) that do not meet the specified QC criteria, creating
#' a new \code{SpatialFeatureExperiment} with the filtered features.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#'
#' @return A new \code{SpatialFeatureExperiment} object with features that pass
#' the QC thresholds.
#'
#' @details This function performs quality control on a
#' \code{SpatialFeatureExperiment} object, filtering out features that do not
#' meet the specified QC criteria. The QC criteria are defined based on the
#' internal \code{.int_filterQCDiscard} function.
#'
#' @seealso \code{\link{.int_filterQCDiscard}} for details on the internal QC
#' filtering function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics feature quality-control
#'
#' @rdname applyQCthresh_loc
#'
#' @aliases applyQCthresh
#'
#' @examples
#' \dontrun{
#' # Apply QC thresholds to a SpatialFeatureExperiment object
#' sfe  # An SFE object for which filters have been defined in colData
#' sfe_filtered <- applyQCthresh_loc(sfe, sample_id = TRUE)
#' }
#'
#' @export
applyQCthresh_loc <- function(sfe, sample_id = TRUE) {
  ## Check arguments
  stopifnot(is(sfe, "SpatialFeatureExperiment"))

  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Create a list of filtered SpatialFeatureExperiment objects
  sfe_list <- lapply(ids, .int_filterQCDiscard, sfe = sfe)

  ## Merge the list of SpatialFeatureExperiment objects into one
  sfeOut <- Reduce(function(x, y) cbind(x, y), sfe_list)

  ## Clean up duplicate entries in the metadata slot
  sfeOut <- .int_cleanMetaData(sfeOut = sfeOut)

  return(sfeOut)
}


#' Apply Quality Control Thresholds to MetaSpatialFeatureExperiment
#'
#' This function applies quality control (QC) thresholds to a
#' \code{MetaSpatialFeatureExperiment} object. It filters out features (genes
#' or other genomic elements) that do not meet the specified QC criteria,
#' creating a new \code{MetaSpatialFeatureExperiment} with the filtered
#' features.
#'
#' @param msfe A SpatialFeatureExperiment or a MetaSpatialFeatureExperiment.
#'
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#'
#' @return A new \code{MetaSpatialFeatureExperiment} object with features that
#' pass the QC thresholds.
#'
#' @details This function performs quality control on a
#' \code{MetaSpatialFeatureExperiment} object, filtering out features that do
#' not meet the specified QC criteria. The QC criteria are defined based on the
#' internal \code{.int_discardFeat} function.
#'
#' @seealso \code{\link{.int_discardFeat}} for details on the internal QC
#' filtering function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics feature quality-control
#'
#' @rdname applyQCthresh_feat
#'
#' @examples
#' \dontrun{# Apply QC thresholds to a MetaSpatialFeatureExperiment object
#' msfe <- # A metaSFE object for which filters have been defined in rowData
#' msfe_filtered <- applyQCthresh_feat(msfe)
#' }
#'
#' @export
applyQCthresh_feat <- function(msfe,
                               sample_id = TRUE) {
  ## Select samples
  ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)

  ## Filter genes
  msfe_int <- lapply(msfe@sfe_data[ids], .int_discardFeat)

  ## If specific samples where modified replace in the metaSFE list
  if (is.character(sample_id)) {
    msfe@sfe_data[names(msfe_int)] <- msfe_int
  } else {
    msfe@sfe_data <- msfe_int
  }

  return(msfe)
}


#' Internal Function: Filter QC Discard Status for a Specific Sample in
#' SpatialFeatureExperiment
#'
#' This internal function filters out features marked for quality control (QC)
#' discard for a specific sample in a \code{SpatialFeatureExperiment} object.
#'
#' @param id A character string specifying the sample identifier.
#' @param sfe A \code{SpatialFeatureExperiment} object.
#'
#' @return A new \code{SpatialFeatureExperiment} object with QC-discarded
#' features removed for the specified sample.
#'
#' @details This internal function filters out features that have been marked
#' for quality control (QC) discard for a specific sample. It returns a
#' modified \code{SpatialFeatureExperiment} object with the QC-discarded
#' features removed.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal spatial transcriptomics quality-control
#'
#' @rdname dot-int_filterQCDiscard
#'
.int_filterQCDiscard <- function(id, sfe) {
  sfe_int <- sfe[, colData(sfe)$sample_id == id]
  keep <- Matrix::which(colData(sfe_int)$qc_discard == FALSE)
  sfe_int <- sfe_int[, keep]
  return(sfe_int)
}


#' Internal Function: Clean Up Duplicate Entries in the Metadata Slot of a
#' SpatialFeatureExperiment
#'
#' This internal function cleans up duplicate entries in the metadata slot of a
#' \code{SpatialFeatureExperiment} object.
#'
#' @param sfeOut A \code{SpatialFeatureExperiment} object.
#'
#' @return A modified \code{SpatialFeatureExperiment} object with duplicate
#' entries removed from the metadata slot.
#'
#' @details This internal function ensures that the metadata slot of the
#' \code{SpatialFeatureExperiment} object does not contain duplicate entries.
#' It makes the names of the metadata unique and keeps the slots with the
#' original names.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal spatial transcriptomics metadata
#'
#' @rdname dot-int_cleanMetaData
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


#' Internal Function: Discard Features in a SpatialFeatureExperiment Based on
#' Quality Control
#'
#' This internal function discards features in a \code{SpatialFeatureExperiment}
#'  object based on the quality control information stored in the feature
#'  metadata (rowData).
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#'
#' @return A modified \code{SpatialFeatureExperiment} object with discarded
#' features removed.
#'
#' @details This internal function discards features in the
#' \code{SpatialFeatureExperiment} object where the quality control information
#' in the feature metadata, specifically the \code{qc_discard} column, is
#' \code{TRUE}. Features marked for discard are removed from the object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal spatial transcriptomics quality control
#'
#' @rdname dot-int_discardFeat
#'
.int_discardFeat <- function(sfe){
  sfe <- sfe[!rowData(sfe)[["qc_discard"]], ]
  return(sfe)
}
