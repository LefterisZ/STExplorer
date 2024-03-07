#' Get Subset of Features from SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment
#'
#' This function retrieves a subset of features based on the specified criteria
#' from either a SpatialFeatureExperiment (SFE) or a
#' MetaSpatialFeatureExperiment (MetaSFE) object.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (TRUE for all IDs in m_sfe).
#' @param subset A regular expression string to be used for subsetting.
#' @param set Either "rowData" or "colData" specifying whether to operate on
#' row or column data.
#' @param col_name The name of the column to be used for subsetting if set to
#' "rowData" or "colData".
#'
#' @return If m_sfe is a SpatialFeatureExperiment (SFE) object, it returns
#' a logical vector indicating the subset of features based on the criteria.
#' If m_sfe is a MetaSpatialFeatureExperiment (MetaSFE) object, it returns
#' a list where each element is a logical vector indicating the subset of
#' samples for each SFE object within the MetaSFE.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # For SFE object
#' sfe <- SpatialFeatureExperiment()
#' getSubset(sfe, sample_id = TRUE, subset = "(^MT-)|(^mt-)", set = "rowData", col_name = "gene_name")
#'
#' # For MetaSFE object
#' metaSFE <- MetaSpatialFeatureExperiment(list(sfe1, sfe2))
#' getSubset(metaSFE, sample_id = TRUE, subset = "(^MT-)|(^mt-)", set = "colData", col_name = "sample_id")
#' }
#'
#' @rdname getSubset
#'
#' @export
getSubset <- function(m_sfe,
                      sample_id,
                      subset,
                      set = c("colData", "rowData"),
                      col_name = NULL) {
  UseMethod("getSubset")
}

#' @export
getSubset.SpatialFeatureExperiment <- function(m_sfe,
                                               sample_id,
                                               subset, set,
                                               col_name = NULL) {
  if (set == "rowData") {
    if (is.null(col_name)) {
      stop("For 'rowData', 'col_name' must be provided.")
    }
    return(grepl(subset, rowData(m_sfe)[[col_name]]))
  } else if (set == "colData") {
    if (is.null(col_name)) {
      stop("For 'colData', 'col_name' must be provided.")
    }
    return(grepl(subset, colData(m_sfe)[[col_name]]))
  } else {
    stop("Invalid 'set' argument. It should be either 'rowData' or 'colData'.")
  }
}


#' @export
getSubset.MetaSpatialFeatureExperiment <- function(m_sfe,
                                                   sample_id,
                                                   subset,
                                                   set,
                                                   col_name = NULL) {
  if (set == "rowData") {
    if (is.null(col_name)) {
      stop("For 'rowData', 'col_name' must be provided.")
    }
    sfe_list <- lapply(m_sfe@sfe_data,
                       .int_msfeGetSubSetRow,
                       col_name = col_name,
                       subset = subset)
    return(sfe_list)
  } else if (set == "colData") {
    if (is.null(col_name)) {
      stop("For 'colData', 'col_name' must be provided.")
    }
    sfe_list <- lapply(m_sfe@sfe_data,
                       .int_msfeGetSubSetCol,
                       col_name = col_name,
                       subset = subset)
    return(sfe_list)
  } else {
    stop("Invalid 'set' argument. It should be either 'rowData' or 'colData'.")
  }
}


# ---------------------------------------------------------------------------- #
#  ### INTERNAL FUNCTIONS ASSOCIATED WITH RETRIEVING A SUBSET OF FEATURES ###
# ---------------------------------------------------------------------------- #
#' INTERNAL: Retrieve Subset of Columns based on Regular Expression from a
#' SpatialFeatureExperiment
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param col_name The name of the column to be used for subsetting.
#' @param subset A regular expression string to be used for subsetting.
#'
#' @return A logical vector indicating the subset of features based on the
#' criteria.
#'
#' @rdname dot-int_msfeGetSubSetRow
#'
.int_msfeGetSubSetRow <- function(sfe, col_name, subset) {
  return(grepl(subset, rowData(sfe)[[col_name]]))
}

#' INTERNAL: Retrieve Subset of Rows based on Regular Expression from a
#' SpatialFeatureExperiment
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param col_name The name of the column to be used for subsetting.
#' @param subset A regular expression string to be used for subsetting.
#'
#' @return A logical vector indicating the subset of features based on the
#' criteria.
#'
#' @rdname dot-int_msfeGetSubSetCol
#'
.int_msfeGetSubSetCol <- function(sfe, col_name, subset) {
  return(grepl(subset, colData(sfe)[[col_name]]))
}
