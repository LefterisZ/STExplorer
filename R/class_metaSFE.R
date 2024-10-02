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
#' @importFrom methods new
#'
#' @examples
#' ## Create an empty MetaSpatialFeatureExperiment object
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
  x@sample_ids[[id]] <- id
  x
}

#' Method to access a SpatialFeatureExperiment by name
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param sample_id Character string. The name of the experiment/sample for
#' the SFE to be accessed. If left to NULL, then the first sfe object in the
#' list is selected.
#'
#' @return The SFE object associated with the provided name.
#'
#' @export
getSFE <- function(x, sample_id = NULL) {
  if (is.null(sample_id)) {
    x@sfe_data[[1]]
  } else {
    x@sfe_data[[sample_id]]
  }
}

#' Method to add an SFE object with multiple samples to a
#' MetaSpatialFeatureExperiment object
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param sfe An instance of the SpatialFeatureExperiment class with multiple
#' samples.
#'
#' @return An updated MetaSpatialFeatureExperiment object with the added SFEs.
#'
#' @importFrom terra unwrap
#' @importFrom SpatialFeatureExperiment SpatRasterImage
#'
#' @export
addMultiSFE <- function(x, sfe) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = TRUE)

  for (i in seq_along(ids)) {
    id <- ids[i]

    ## If the multi-sample SFE object was saved and loaded again into the
    ## environment, the images are going to be of class "PackedRasterImage". If
    ## this is the case, we need to "unwrap" the packaged images.
    is_Packed <- inherits(imgData(sfe)$data[[i]], "PackedRasterImage")
    if (is_Packed) {
      img <- SpatRasterImage(unwrap(imgData(sfe)$data[[i]]))
      sfe@int_metadata$imgData$data[[i]] <- img
    }

    ## Subset per sample
    sfe_int <- sfe[, colData(sfe)$sample_id == id]

    ## Add each subset to the MetaSpatialFeatureExperiment object
    x <- addSFE(x, sfe_int, sample_id = id)
  }

  x
}

#' Method to retrieve multiple SFE objects as a single SFE object from within
#' an MSFE object
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param sample_ids Character string. The sample IDs from the SFE objects to
#' be combined and retrieved as one SFE object
#'
#' @details
#' During merging the SFE objects, the common columns from the rowData are
#' removed to allow merging. However, the unique rowData columns from each SFE
#' object are kept. Additionally, the "gene_name", "id", "mean", "detected",
#' and "total" columns are generated again. Only this time, the "mean",
#' "detected", and "total" values are calculated based on the gene counts in
#' ALL samples and NOT in a per sample way.
#'
#' The colData are easier to merge since there we are performing an rbind
#' instead of the cbind we have to do in the rowData. As a result all columns
#' are kept.
#'
#' @return An SFE object with multiple samples inside.
#'
#' @importFrom SpatialFeatureExperiment cbind
#'
#' @export
getMultiSFE <- function(x, sample_ids) {
  ## Check provided sample IDs are correct
  .int_checkSampleNames(x, sample_ids)

  ## Check that all sfe objects have the counts assay as 'dgCMatrix'.
  ## If NOT, then transform to 'dgCMatrix'.
  for (id in sample_ids) {
    sfe <- getSFE(x, id)
    mtx_nm <- names(SummarizedExperiment::assays(sfe))

    for (nm in mtx_nm) {
      if (!is(assay(sfe, nm), "dgCMatrix")) {
        sfe <- .int_convertToDgCMatrix(sfe, nm)
      }
    }

    ## Place the updated sfe object back into the msfe
    x@sfe_data[[id]] <- sfe
  }

  ## Get some data to place back into the rowData
  gene_name <- rowData(x@sfe_data[[1]])[["gene_name"]]

  ## Extract all rowData from sfe objects
  all_row_data <- lapply(x@sfe_data, rowData)

  ## Determine common columns
  # Assuming all sfe objects have the same columns for simplicity
  # If not, we'll need to find the intersection of column names
  common_cols <- Reduce(intersect, lapply(all_row_data, colnames))

  ## Remove common columns from each sfe object's rowData
  x@sfe_data <- lapply(x@sfe_data, function(sfe) {
    # Only keep columns that are not common
    rowData(sfe) <- rowData(sfe)[, !colnames(rowData(sfe)) %in% common_cols, drop = FALSE]
    return(sfe)
  })

  ## Extract SFE objects based on sample_ids
  sfe_list <- lapply(sample_ids, function(id) getSFE(x, id))

  ## Use Reduce with cbind to combine SFE objects
  sfe_out <- Reduce(SpatialFeatureExperiment::cbind, sfe_list)

  ## Export the unique columns to put them back after we add the gene-QC metrics
  ## Otherwise it interferes with the columns added by `scater::addPerFeatureQC`
  rowdata <- rowData(sfe_out)

  ## Empty the rowData from the sfe_out
  rowData(sfe_out) <- NULL

  ## Place columns back into rowData
  rowData(sfe_out)[["symbol"]] <- gene_name

  ## Add basic gene-QC metrics --> but now calculated for all samples as one
  sfe_out <- addPerGeneQC(sfe_out,
                          sample_id = id,
                          assay = "counts",
                          version = NULL,
                          mirror = NULL,
                          add = "none")

  ## Place back the unique columns in the rowData
  rowData(sfe_out) <- merge(rowData(sfe_out),
                            rowdata,
                            by = "row.names",
                            sort = FALSE) %>%
    data.frame() %>%
    column_to_rownames(var = "Row.names") %>%
    DataFrame()

  return(sfe_out)
}

#' Method to retrieve sample IDs from MetaSpatialFeatureExperiment
#'
#' This method retrieves the sample IDs from a MetaSpatialFeatureExperiment
#' object.
#'
#' The sample IDs are stored in the @sample_ids slot of the
#' MetaSpatialFeatureExperiment.
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
getSampleIDs <- function(x) {
  unlist(x@sample_ids)
}


# ---------------------------------------------------------------------------- #
#  ########### INTERNAL FUNCTIONS ASSOCIATED WITH THE MSFE CLASS ############
# ---------------------------------------------------------------------------- #
#' Internal: Check if sample names are correct
#'
#' An internal function to check that the sample names provided for extraction
#' are correct.
#'
#' @param x An instance of the MetaSpatialFeatureExperiment class.
#' @param sample_ids Character string. The sample IDs to be checked.
#'
#' @rdname dot-int_checkSampleNames
#'
#' @return Logical vector indicating whether each sample ID is correct.
.int_checkSampleNames <- function(x, sample_ids) {
  test_names <- sample_ids %in% getSampleIDs(x)
  if (!all(test_names)) {
    no_match <- sample_ids[!test_names]
    stop("Some sample names do not match the names in the MSFE object ",
         "provided in the `x` argument.\n",
         "Please check the names you provided in the `sample_ids` argument.\n",
         "The non-matching names are: \n",
         paste(no_match, collapse = ", "))
  }
}

#' Internal: Convert assay matrix to dgCMatrix
#'
#' An internal function to convert any assay matrix class to dgCMatrix class.
#'
#' @param sfe An instance of the SummarizedExperiment class.
#' @param nm Name of the assay matrix.
#'
#' @rdname dot-int_convertToDgCMatrix
#'
#' @return SummarizedExperiment with the specified assay matrix converted to
#' dgCMatrix.
.int_convertToDgCMatrix <- function(sfe, nm) {
  assay(sfe, nm) <- as(assay(sfe, nm), "dgCMatrix")
  return(sfe)
}

#' Internal Function: Get SpatialFeatureExperiment Sample IDs
#'
#' This internal function retrieves the sample IDs based on the provided
#' criteria.
#'
#' @param msfe A SpatialFeatureExperiment object.
#' @param sample_id Either a logical vector indicating the samples to include
#' (if TRUE, all samples are included), NULL (if NULL selects the first sample),
#' or a character vector specifying the sample IDs to include.
#'
#' @return A character vector of selected sample IDs.
#'
.int_getMSFEsmplID <- function(msfe, sample_id) {
  ## Select samples
  if (isTRUE(sample_id)) {
    ids <- names(msfe@sfe_data)
  } else if (is.character(sample_id)) {
    ids <- sample_id
  } else if (is.null(sample_id)) {
    ids <- names(msfe@sfe_data)[1]
  }

  return(ids)
}


#' Internal Function: Get SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment
#'
#' This internal function checks the class type of the input object and returns
#' either the original SpatialFeatureExperiment (SFE) object or, in the case of
#' a MetaSpatialFeatureExperiment (MetaSFE), retrieves a subset based on the
#' provided sample ID.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant for MetaSpatialFeatureExperiment).
#'
#' @return An object of class SpatialFeatureExperiment. If the input is already
#' an SFE, it returns the original object. If the input is a MetaSFE, it
#' returns a subset based on the provided sample IDs.
#'
#' @rdname dot-int_sfeORmsfe
#'
.int_sfeORmsfe <- function(m_sfe, sample_id) {
  SFE <- is(m_sfe, "SpatialFeatureExperiment")
  metaSFE <- is(m_sfe, "MetaSpatialFeatureExperiment")
  if (SFE) {
    sfe <- m_sfe
  }
  if (metaSFE) {
    sfe <- getSFE(m_sfe, sample_id)
  }

  return(sfe)
}


#' Internal: Check and Output MSFE or SFE Object
#'
#' An internal function to check the class type of the input object
#' (\code{m_sfe}) and return either the original SpatialFeatureExperiment (SFE)
#' object or, in the case of a MetaSpatialFeatureExperiment (MetaSFE), add a
#' subset based on the provided sample IDs (\code{sample_id}).
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sfe An object of class SpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample ID for the sample
#' we used.
#'
#' @return An object of class SpatialFeatureExperiment. If the input is already
#' an SFE, it returns the original object. If the input is a MetaSFE, it
#' returns a subset based on the provided sample IDs.
#'
#' @seealso \code{\link{addSFE}}
#'
#' @rdname dot-int_checkAndOutput
.int_checkAndOutput <- function(m_sfe, sfe, sample_id) {
  # Check if m_sfe is a MetaSpatialFeatureExperiment object
  MSFE <- is(m_sfe, "MetaSpatialFeatureExperiment")

  # If MSFE, add SFE using the addSFE function
  if (MSFE) {
    out <- addSFE(x = m_sfe, sfe = sfe, sample_id = sample_id)
  } else {
    out <- sfe
  }

  # Return the resulting object
  return(out)
}

