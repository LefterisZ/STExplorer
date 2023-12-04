#' Set Quality Control Threshold for Library Size
#'
#' This function sets quality control thresholds for library size in a
#' SpatialFeatureExperiment object.
#'
#' @param sfe A SpatialFeatureExperiment object containing count data.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param min_t Numeric, the minimum library size threshold. If NA, no lower
#' threshold is applied. Default is NA.
#' @param max_t Numeric, the maximum library size threshold. If NA, no upper
#' threshold is applied. Default is NA.
#'
#' @return A SpatialFeatureExperiment object with a new column 'qc_lib_size'
#' indicating whether the library size passes the thresholds.
#'
#' @details This function filters out locations based on the specified library
#' size thresholds. It also prints out a summary table showing the number of filtered
#' locations per sample.
#'
#' @seealso \code{\link{.int_getSmplIDs}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords data manipulation
#'
#' @rdname setQCthresh_LibSize
#'
#' @examples
#' # 'sfe' is a SpatialFeatureExoeriment object
#'
#' # Both filters
#' sfe <- setQCthresh_LibSize(sfe, sample_id = "JBO019", 500, 35000)
#'
#' # Upper only
#' sfe <- setQCthresh_LibSize(sfe, sample_id = "JBO019", NA, 35000)
#'
#' # Lower only
#' sfe <- setQCthresh_LibSize(sfe, sample_id = "JBO019", 500, NA)
#'
#' @export
setQCthresh_LibSize <- function(sfe, sample_id = TRUE, min_t, max_t) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Select threshold
  if (is.na(min_t)) {
    thresh_min <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_min <- colData(sfe)$sum < min_t
  }
  if (is.na(max_t)) {
    thresh_max <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_max <- colData(sfe)$sum > max_t
  }
  qc_lib_size <- (thresh_min | thresh_max) & (colData(sfe)$sample_id %in% ids)

  ## Check how many spots are filtered out
  cat("Number of locations filtered out:\n")
  print(table(colData(sfe)$sample_id, qc_lib_size))

  ## Add threshold in colData
  colData(sfe)$qc_lib_size <- qc_lib_size

  return(sfe)
}


#' Set QC Threshold for Genes Expressed in each location
#'
#' This function sets quality control (QC) thresholds for the number of
#' expressed genes in each location.
#'
#' @param sfe A SpatialFeatureExperiment object containing gene expression data.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param min_t Numeric, the minimum threshold for detection. Locations with
#' number of expressed genes below this threshold will be flagged. Use NA for
#' no minimum threshold.
#' @param max_t Numeric, the maximum threshold for detection. Locations with
#' number of expressed genes above this threshold will be flagged. Use NA for
#' no maximum threshold.
#'
#' @return A modified SpatialFeatureExperiment object with an added column
#' indicating whether each location passes the QC thresholds.
#'
#' @details The function sets QC thresholds for the number of expressed genes
#' in each location based on the specified minimum and maximum detection values.
#'  Locations failing to meet the thresholds are flagged, and a summary table
#'  of filtered locations is printed.
#'
#' @seealso \code{\link{.int_getSmplIDs}}, \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, genes expressed, SpatialFeatureExoeriment
#'
#' @rdname setQCthresh_GenesExpr
#'
#' @aliases setQCthresh_GenesExpr
#'
#' @examples
#' # 'sfe' is a SpatialFeatureExoeriment object
#'
#' # Set QC thresholds for number of expressed genes in each location
#' sfe <- setQCthresh_GenesExpr(sfe, sample_id = "JBO019", min_t = 100)
#'
#' @export
setQCthresh_GenesExpr <- function(sfe, sample_id = TRUE, min_t, max_t) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Select threshold
  if (is.na(min_t)) {
    thresh_min <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_min <- colData(sfe)$detected < min_t
  }
  if (is.na(max_t)) {
    thresh_max <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_max <- colData(sfe)$detected > max_t
  }
  qc_detected <- (thresh_min | thresh_max) & (colData(sfe)$sample_id %in% ids)

  ## Check how many spots are filtered out
  cat("Number of locations filtered out:\n")
  print(table(colData(sfe)$sample_id, qc_detected))

  ## Add threshold in colData
  colData(sfe)$qc_detected <- qc_detected

  return(sfe)
}


#' Set QC Threshold for Mitochondrial Content
#'
#' This function sets quality control (QC) thresholds for mitochondrial content
#' based on the percentage of mitochondrial gene expression form the total of
#' expressed genes in a location.
#'
#' @param sfe A SpatialFeatureExperiment object containing mitochondrial
#' content data (a colData column named subset_mito_percent)
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param min_t Numeric, the minimum threshold for mitochondrial content
#' percentage. Samples with percentages below this threshold will be flagged.
#' Use NA for no minimum threshold.
#' @param max_t Numeric, the maximum threshold for mitochondrial content
#' percentage. Samples with percentages above this threshold will be flagged.
#' Use NA for no maximum threshold.
#'
#' @return A modified SpatialFeatureExperiment object with an added column in
#' the colData indicating whether each location passes the QC thresholds for
#' mitochondrial content.
#'
#' @details The function sets QC thresholds for mitochondrial content based on
#' the specified minimum and maximum percentage values. Locations failing to
#' meet the thresholds are flagged, and a summary table of filtered locations
#' is printed.
#'
#' @seealso \code{\link{.int_getSmplIDs}}, \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, mitochondrial content, SpatialFeatureExperiment
#'
#' @rdname setQCthresh_Mito
#'
#' @aliases setQCthresh_Mito
#'
#' @examples
#' # 'sfe' is a SpatialFeatureExoeriment object
#'
#' # Set QC thresholds for mitochondrial content percentages above 30
#' sfe <- setQCthresh_Mito(sfe, sample_id = TRUE, max_t = 30)
#'
#' @export
setQCthresh_Mito <- function(sfe, sample_id = TRUE, min_t, max_t) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Select threshold
  if (is.na(min_t)) {
    thresh_min <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_min <- colData(sfe)$subsets_mito_percent < min_t
  }
  if (is.na(max_t)) {
    thresh_max <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_max <- colData(sfe)$subsets_mito_percent > max_t
  }
  qc_mito <- (thresh_min | thresh_max) & (colData(sfe)$sample_id %in% ids)

  ## Check how many spots are filtered out
  cat("Number of locations filtered out:\n")
  print(table(colData(sfe)$sample_id, qc_mito))

  ## Add threshold in colData
  colData(sfe)$qc_mito <- qc_mito

  return(sfe)
}


#' Set QC Threshold for Cell Counts
#'
#' This function sets quality control (QC) thresholds for cell counts based on
#' the values in the cellCount column.
#'
#' @param sfe A SpatialFeatureExoeriment object containing cell count data.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param min_t Numeric, the minimum threshold for cell counts. Locations with
#' counts below this threshold will be flagged. Use NA for no minimum threshold.
#' @param max_t Numeric, the maximum threshold for cell counts. Locations with
#' counts above this threshold will be flagged. Use NA for no maximum threshold.
#'
#' @return A modified SpatialFeatureExoeriment object with an added column in
#' colData indicating whether each location passes the QC thresholds for cell
#' counts.
#'
#' @details The function sets QC thresholds for cell counts based on the
#' specified minimum and maximum values. Locations failing to meet the
#' thresholds are flagged, and a summary table of filtered locations is printed.
#'
#' @seealso \code{\link{.int_getSmplIDs}}, \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, cell counts, SpatialFeatureExoeriment
#'
#' @rdname setQCthresh_CellCount
#'
#' @aliases setQCthresh_CellCount
#'
#' @examples
#' # 'sfe' is a SpatialFeatureExoeriment object
#'
#' # Set QC thresholds with sample IDs for cell counts above 15
#' sfe <- setQCthresh_CellCount(sfe, sample_id = TRUE, max_t = 15)
#'
#' @export
setQCthresh_CellCount <- function(sfe, sample_id = TRUE, min_t, max_t) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Select threshold
  if (is.na(min_t)) {
    thresh_min <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_min <- colData(sfe)$cellCount < min_t
  }
  if (is.na(max_t)) {
    thresh_max <- rep(FALSE, n = ncol(sfe))
  } else {
    thresh_max <- colData(sfe)$cellCount > max_t
  }
  qc_cellCount <- (thresh_min | thresh_max) & (colData(sfe)$sample_id %in% ids)

  ## Check how many spots are filtered out
  cat("Number of locations filtered out:\n")
  print(table(colData(sfe)$sample_id, qc_cellCount))

  ## Add threshold in colData
  colData(sfe)$qc_cellCount <- qc_cellCount

  return(sfe)
}


#' Set QC Threshold for Missing Annotations
#'
#' This function sets quality control (QC) thresholds for missing annotations
#' based on the presence of NAs in the annotation column. This QC function is
#' mostly useful when there are relatively few parts of the tissue sample that
#' couldn't be annotated by an expert or another approach.
#'
#' @param sfe A SpatialFeatureExoeriment object containing annotation data.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#'
#' @return A modified SpatialFeatureExoeriment object with an added column  in
#' colData indicating whether each location has missing annotations.
#'
#' @details The function sets QC thresholds for missing annotations based on
#' the presence of NAs in the annotation column in colData. Locations with
#' missing annotations are flagged, and a summary table of filtered locations
#' is printed. his QC function is
#' mostly useful when there are relatively few parts of the tissue sample that
#' couldn't be annotated by an expert or another approach.
#'
#' @seealso \code{\link{.int_getSmplIDs}}, \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, missing annotations, SpatialFeatureExoeriment
#'
#' @rdname setQCthresh_NAs
#'
#' @aliases setQCthresh_NAs
#'
#' @examples
#' # 'sfe' is a SpatialFeatureExoeriment object
#'
#' # Set QC thresholds with for missing annotations
#' sfe <- setQCthresh_NAs(sfe, sample_id = TRUE)
#'
#' @export
setQCthresh_NAs <- function(sfe, sample_id = TRUE) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  is_na <- is.na(colData(sfe)$annotation)

  qc_NA_spots <- is_na & (colData(sfe)$sample_id %in% ids)

  ## Check how many spots are filtered out
  cat("Number of locations filtered out:\n")
  print(table(colData(sfe)$sample_id, qc_NA_spots))

  ## Add threshold in colData
  colData(sfe)$qc_NA_spots <- qc_NA_spots

  return(sfe)
}

#' Select Locations to Discard
#'
#' This function selects locations to discard based on specified QC metrics.
#'
#' @param sfe A SpatialFeatureExoeriment object containing QC metrics.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param filters Either TRUE to select all calculated QC metrics or a
#' character vector with the column names of the QC metrics to consider for
#' filtering locations. QC metric columns begin with "qc_".
#'
#' @return A modified SpatialFeatureExoeriment object with an added column indicating whether each location passes the specified QC thresholds.
#'
#' @details The function sets QC thresholds for discarding locations based on the specified QC metrics. Locations failing to meet the thresholds are flagged, and a summary table of filtered locations is printed. Users can choose to filter based on all available QC metrics or specify a subset using the `filters` argument.
#'
#' @seealso \code{\link{.int_getSmplIDs}}, \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, location filtering, SpatialFeatureExoeriment
#'
#' @rdname setQCtoDiscard_loc
#'
#' @aliases setQCtoDiscard_loc
#'
#' @examples
#' # Set QC thresholds with sample IDs for discarding locations based on all available QC metrics
#' se <- setQCtoDiscard_loc(sfe, sample_id = TRUE, filters = TRUE)
#'
#' @export
setQCtoDiscard_loc <- function(sfe, sample_id = TRUE, filters = TRUE) {
  ## Get sample IDs
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Select QC metrics to use
  if (isTRUE(filters)) {
    # Specify columns
    selected_cols <- grepl("qc_", colnames(colData(sfe)))

  } else if (is.character(filters)) {
    # Concatenate the elements with "|" separator
    pattern <- paste(filters, collapse = "|")
    # Specify columns
    selected_cols <- grepl(pattern, colnames(colData(sfe)))
    # Check all columns exist
    exist <- sum(selected_cols) == length(filters)
    if (exist) {
      message("Using user-defined QC columns: ",
              paste(filters, collapse = " "))
    } else {
      stop("\nNot all QC columns defined by the user are present in the ",
           "colData(sfe) column names.\n",
           "Column names provided: ", paste(filters, collapse = " "),
           "\nExisting QC columns: ",
           paste(colnames(colData(sfe))[grepl("qc_", colnames(colData(sfe)))],
                 collapse = " "))
    }

  } else {
    stop("\n- The filters argument should be either TRUE or a character \n",
         "vector containing the names of the QC columns you require.\n",
         "- The QC columns start with 'qc_'.\n",
         "- You can find them typing: colnames(colData(sfe)) \n")
  }

  ## Create a vector where each element corresponds to a row in the data frame
  qc_discard <- apply(colData(sfe)[selected_cols], 1, function(row) any(row))
  names(qc_discard) <- NULL

  ## Make sure the result has only logical (TRUE/FALSE) values
  qc_discard <- as.logical(qc_discard)

  ## Check the number of discarded spots for each metric
  cat("Number of discarded locations for each QC metric:\n")
  print(apply(colData(sfe)[selected_cols], 2, sum))

  ## Check how many spots are filtered out
  cat("Number of locations filtered out:\n")
  print(table(colData(sfe)$sample_id, qc_discard))

  ## Add threshold in colData
  colData(sfe)$qc_discard <- qc_discard

  return(sfe)
}


#' Set QC Thresholds for Zero Expression
#'
#' These functions set quality control (QC) thresholds for gene expression data
#' based on zero expression values.
#'
#' @param msfe A SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#'
#' @return A SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' object with added columns in the rowData indicating whether each gene passes
#' the specified QC thresholds.
#'
#' @details This function, identifies genes without expression in each sample.
#'
#' @seealso \code{\link{.int_isZero}}, \code{\link{.int_logLowMean}},
#' \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, gene expression, SpatialFeatureExoeriment
#'
#' @rdname setQCthresh_ZeroExpr
#'
#' @aliases setQCthresh_ZeroExpr
#'
#' @examples
#' # 'sfe' is a SpatialFeatureExoeriment object
#'
#' # Set QC thresholds for zero expression
#' msfe <- setQCthresh_ZeroExpr(msfe, sample_id = TRUE)
#'
#' @export
setQCthresh_ZeroExpr <- function(msfe, sample_id = TRUE) {
  ## Find genes without expression in each sample
  msfe_int <- lapply(msfe, .int_isZero)

  return(msfe_int)
}


#' Set QC Thresholds for Low Logarithmic Mean in Gene Expression
#'
#' This function sets quality control (QC) thresholds for low logarithmic mean
#' values in gene expression data. It is used after log transforming the gene
#' count matrix
#'
#' @param msfe A SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param threshold Numeric, the threshold for the logarithmic mean value.
#' Genes with mean values below this threshold will be flagged.
#'
#' @return A SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' object with added column in rowData indicating whether each gene passes the
#' specified QC threshold for low logarithmic mean values.
#'
#' @details The function sets QC thresholds for low logarithmic mean values in
#' gene expression data for each sample in the provided object. Genes failing
#' to meet the thresholds are flagged.
#'
#' @seealso \code{\link{.int_logLowMean}}, \code{\link[SpatialFeatureExperiment]{crowData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, gene expression, SpatialFeatureExoeriment
#'
#' @rdname setQCthresh_LowLogMean
#'
#' @aliases setQCthresh_LowLogMean
#'
#' @examples
#' # 'msfe' is a SpatialFeatureExoeriment or a MetaSpatialFeatureExoeriment,
#' # object
#'
#' # Set QC thresholds for low logarithmic mean values of gene expression.
#' msfe <- setQCthresh_LowLogMean(msfe, threshold = 1, sample_id = TRUE)
#'
#' @export
setQCthresh_LowLogMean <- function(msfe,
                                   threshold = 1,
                                   sample_id = TRUE) {
  ## Get sample IDs
  ids <- names(msfe)

  ## Initialise list
  msfe_int <- list()

  ## Perform calculations
  for (id in ids) {
    msfe_int[[id]] <- .int_logLowMean(sfe = msfe[[id]],
                                      threshold = threshold)
  }

  return(msfe_int)
}


#' Set QC Thresholds using Custom Metrics
#'
#' This function sets quality control (QC) thresholds using custom metrics for
#' gene expression data or for location data.
#'
#' @param msfe A SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param MARGIN Numeric, indicating whether to apply the custom metric along
#' rows (1) or columns (2). If rows is selected then the filters are applied to
#' genes and the flagged genes are stored in a column in rowData. If columns is
#' selected then the filters are applied to locations and the flagged locations
#' are stored in a column in colData.
#' @param qcMetric A vector of TRUE and FALSE for genes or locations that are
#' going to be discarded. You need to have prepared this vector in advance. The
#' vector needs to be of length equal to nrow(rowData(sfe)) if MARGIN == 1, or
#' equal to nrow(colData(sfe)) if MARGIN == 2.
#'
#' @return A SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' object with added columns indicating whether each gene or location passes
#' the specified QC thresholds using custom metrics.
#'
#' @details The function sets QC thresholds for gene expression data using
#' custom metrics provided in the `qcMetric` vector. The user can specify
#' whether to apply the custom metric along rows (genes) or columns (samples)
#' using the `MARGIN` argument. Genes or locations failing to meet the
#' thresholds are flagged.
#'
#' @seealso \code{\link{.int_custom_1}}, \code{\link{.int_custom_2}},
#' \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, gene expression, SummarizedExperiment
#'
#' @rdname setQCthresh_custom
#'
#' @aliases setQCthresh_custom
#'
#' @examples
#' # 'msfe' is a SpatialFeatureExoeriment or a MetaSpatialFeatureExoeriment,
#' # object
#'
#' # Set QC thresholds using custom metrics for specific samples
#' qcMetric <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
#' msfe <- setQCthresh_custom(msfe,
#'                            sample_id = c("JBO019"),
#'                            MARGIN = 1,
#'                            qcMetric)
#'
#' @export
setQCthresh_custom <- function(msfe, sample_id = TRUE, MARGIN, qcMetric) {
  ## Select samples
  if (isTRUE(sample_id)) {
    ids <- names(msfe)
  } else if (is.character(sample_id)) {
    ids <- sample_id
  }

  ## Get QC name
  qc_name <- paste0("qc_", deparse(substitute(qcMetric)))

  ## Add QC metric
  if (MARGIN == 1) {
    msfe_int <- lapply(msfe[ids], .int_custom_1, qcMetric, qc_name)
  } else if (MARGIN == 2) {
    msfe_int <- lapply(msfe[ids], .int_custom_2, qcMetric, qc_name)
  }

  ## If specific samples where modified replace in the metaSFE list
  if (is.character(sample_id)) {
    msfe[names(msfe_int)] <- msfe_int
  } else {
    msfe <- msfe_int
  }

  return(msfe)
}


#' Set QC Thresholds for Discarding Features
#'
#' This function marks the features (genes) flagged for discarding based on
#' specified QC metrics.
#'
#' @param msfe A SpatialFeatureExoeriment or a MetaSpatialFeatureExoeriment.
#' @param filters Either TRUE to select all calculated QC metrics or a
#' character vector with the column names of the QC metrics to consider for
#' filtering features. QC metric columns begin with "qc_".
#'
#' @return A SpatialFeatureExoeriment or a MetaSpatialFeatureExoeriment object
#' with features flagged for discard based on specified QC thresholds.
#'
#' @details The function marks the features (genes) flagged for discarding
#' based on the specified QC metrics. A feature failing to meet any of the
#' thresholds applied earlier is flagged.
#'
#' A summary table of filtered features is printed. Users can choose to filter
#' based on all available QC metrics or specify a subset using the `filters`
#' argument.
#'
#' @seealso \code{\link{.int_featToDiscard}}, \code{\link{.int_summaryTable}},
#' \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC, feature filtering, SpatialFeatureExoeriment
#'
#' @rdname setQCtoDiscard_feat
#'
#' @aliases setQCtoDiscard_feat
#'
#' @examples
#' # 'msfe' is a SpatialFeatureExoeriment or a MetaSpatialFeatureExoeriment,
#' # object
#'
#' # Mark the features (genes) flagged for discarding
#' msfe <- setQCtoDiscard_feat(msfe, filters = TRUE)
#'
#' @export
setQCtoDiscard_feat <- function(msfe, filters = TRUE) {
  ## Mark features to discard
  msfe_int <- lapply(msfe, .int_featToDiscard, filters = filters)

  ## Get a summary table
  tbl_list <- lapply(msfe_int, .int_summaryTable, filters = filters)

  ## Check the number of discarded features for each metric per sample
  cat("Number of locations filtered out:\n")
  print(rlist::list.rbind(tbl_list))

  return(msfe_int)
}


#' Select QC Metrics for Filtering
#'
#' This function selects QC metrics for filtering based on specified criteria.
#'
#' @param sfe A SpatialFeatureExperiment object containing spatial
#' transcriptomics experiment data.
#' @param filters TRUE selects all calculated QC metrics. Alternatively,
#' provide a character vector with the column names of the QC metrics you want
#' to consider for filtering. The QC metric columns begin with "qc_".
#'
#' @return A character vector representing the selected QC metrics for
#' filtering.
#'
#' @details The function provides flexibility in selecting QC metrics for
#' filtering. If `filters` is set to TRUE, all QC metrics are selected. If it
#' is a character vector, the function validates and selects the specified QC
#' metrics. The selected columns are then returned for further processing.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC metrics, filtering, selection
#'
#' @rdname dot-int_selectFilters
#'
.int_selectFilters <- function(sfe, filters) {
  ## Select QC metrics to use
  if (isTRUE(filters)) {
    # Specify columns
    col_ns <- colnames(rowData(sfe))
    selected_cols <- grep("^qc_", col_ns, value = TRUE)

  } else if (is.character(filters)) {
    # Concatenate the elements with "|" separator
    pattern <- paste(filters, collapse = "|")
    # Specify columns
    selected_cols <- grepl(pattern, colnames(rowData(sfe)))
    # Check all columns exist
    exist <- sum(selected_cols) == length(filters)
    if (exist) {
      message("Using user-defined QC columns: ",
              paste(filters, collapse = " "))
    } else {
      stop("\nNot all QC columns defined by the user are present in the ",
           "colData(sfe) column names.\n",
           "Column names provided: ", paste(filters, collapse = " "),
           "\nExisting QC columns: ",
           paste(colnames(rowData(sfe))[grepl("qc_", colnames(rowData(sfe)))],
                 collapse = " "))
    }

  } else {
    stop("\n- The filters argument should be either TRUE or a character \n",
         "  vector containing the names of the QC columns you require.\n",
         "- The QC columns start with 'qc_'.\n",
         "- You can find them typing: colnames(rowData(sfe)) \n")
  }

  return(selected_cols)
}


#' Feature Grouping: Mark Features to Discard
#'
#' This function marks features (genes) in a spatial transcriptomics experiment
#' to be discarded based on specified QC metrics.
#'
#' @param sfe A SpatialFeatureExperiment containing spatial
#' transcriptomics experiment data.
#' @param filters TRUE selects all calculated QC metrics. Alternatively,
#' provide a character vector with the column names of the QC metrics you want
#' to consider for filtering. The QC metric columns begin with "qc_".
#'
#' @return An updated SpatialFeatureExperiment experiment object with a new
#' column "qc_discard" in the rowData, indicating whether each feature should
#' be discarded based on the selected QC metrics.
#'
#' @details The function uses the internal function `.int_selectFilters` to
#' obtain the selected QC metrics. It then applies a logical operation to mark
#' features to be discarded based on the presence of TRUE values in at least
#' one of the selected QC metrics. The result is added as a new column
#' "qc_discard" in the rowData of the SpatialFeatureExperiment object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords feature grouping, QC metrics, discard, SpatialFeatureExperiment
#'
#' @rdname dot-int_featToDiscard
#'
.int_featToDiscard <- function(sfe, filters) {
  ## Select QC metrics to use
  selected_cols <- .int_selectFilters(sfe = sfe, filters = filters)

  ## Create a vector where each element corresponds to a row in the data frame
  ##  - find all the genes (rows) with TRUE in at least one of the QC metrics
  qc_discard <- apply(rowData(sfe)[selected_cols], 1, function(row) any(row))
  names(qc_discard) <- NULL

  ## Make sure the result has only logical (TRUE/FALSE) values
  qc_discard <- as.logical(qc_discard)

  ## Add the discard column into the rowData
  rowData(sfe)[["qc_discard"]] <- qc_discard

  return(sfe)
}


#' Generate Summary Table for QC Metrics
#'
#' This function generates a summary table for selected QC metrics in a spatial
#' transcriptomics experiment.
#'
#' @param sfe A SpatialFeatureExperiment containing spatial
#' transcriptomics experiment data.
#' @param filters TRUE selects all calculated QC metrics. Alternatively,
#' provide a character vector with the column names of the QC metrics you want
#' to consider for filtering. The QC metric columns begin with "qc_".
#'
#' @return A data frame representing a summary table for the selected QC
#' metrics, showing the sum of each metric across all samples.
#'
#' @details The function uses the internal function `.int_selectFilters` to
#' obtain the selected QC metrics. It then creates a matrix with one row and
#' columns corresponding to the selected metrics. The function calculates the
#' sum of each metric across all samples and populates the summary table
#' accordingly.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords QC metrics, summary table, SpatialFeatureExperiment
#'
#' @rdname dot-int_summaryTable
#'
.int_summaryTable <- function(sfe, filters) {
  selected_cols <- .int_selectFilters(sfe = sfe, filters = filters)
  n_samples <- 1
  n_cols <- length(selected_cols)
  df <- matrix(data = NA, nrow = n_samples, ncol = n_cols) %>%
    as.data.frame()
  colnames(df) <- selected_cols
  df[1, ] <- apply(rowData(sfe)[selected_cols], 2, sum)

  return(df)
}


#' Identify Genes with Zero Expression
#'
#' This function identifies genes with zero expression in a spatial
#' transcriptomics experiment
#'
#' @param sfe A SpatialFeatureExperiment containing spatial
#' transcriptomics experiment data.
#'
#' @return An updated SpatialFeatureExperiment object with a new column
#' "qc_isZero" in the rowData, indicating whether each gene has zero expression.
#'
#' @details The function checks if the total expression for each gene is equal
#' to zero and adds a new column "qc_isZero" in the rowData of the
#' SpatialFeatureExperiment to mark genes with zero expression.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords zero expression, genes, SpatialFeatureExperiment
#'
#' @rdname dot-int_isZero
#'
.int_isZero <- function(sfe) {
  is_zero <- rowData(sfe)$total == 0
  rowData(sfe)$qc_isZero <- is_zero
  sfe
}


#' Identify Genes with Low Logarithmic Mean Expression
#'
#' This function identifies genes with low logarithmic mean expression in a
#' spatial transcriptomics experiment
#'
#' @param sfe A SpatialFeatureExperiment containing spatial
#' transcriptomics experiment data.
#' @param threshold The threshold value for determining low logarithmic mean
#' expression.
#'
#' @return An updated SpatialFeatureExperiment object with a new column
#' "qc_isLogLow" in the rowData, indicating whether each gene has a logarithmic
#' mean expression below the specified threshold.
#'
#' @details The function calculates the logarithmic mean expression for each
#' gene and checks if it is below the specified threshold. It then adds a new
#' column "qc_isLogLow" in the rowData of the single-cell experiment to mark
#' genes with low logarithmic mean expression.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords logarithmic mean expression, genes, SpatialFeatureExperiment
#'
#' @rdname dot-int_logLowMean
#'
.int_logLowMean <- function(sfe, threshold) {
  colnameLogM <- "s_logMean"
  colnameLogL <- "qc_isLogLow"
  is_logLow <- rowData(sfe)[[colnameLogM]] <= threshold
  rowData(sfe)[[colnameLogL]] <- is_logLow
  return(sfe)
}


#' Apply Custom QC Metric to Rows in a SpatialFeatureExperiment object
#'
#' This function adds a custom QC metric to rows (genes) in a
#' spatial transcriptomics experiment.
#'
#' @param sfe A SpatialFeatureExperiment containing spatial
#' transcriptomics experiment data.
#' @param qcMetric A vector of TRUE and FALSE for rows (genes)
#' that are flagged for discarding.
#' @param qc_name The name of the custom QC metric to be added to the rowData.
#'
#' @return An updated SpatialFeatureExperiment object with the custom QC metric
#' added to the rowData.
#'
#' @details The function takes a vector `qcMetric` containing TRUE and FALSE
#' values, representing rows that should be flagged for discarding. It adds a
#' new column with the specified name `qc_name` in the rowData of the
#' SpatialFeatureExperiment object, marking rows based on the custom QC metric.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords custom QC metric, SpatialFeatureExperiment
#'
#' @rdname dot-int_custom_1
#'
.int_custom_1 <- function(sfe, qcMetric, qc_name) {
  rowData(sfe)[[qc_name]] <- qcMetric

  return(sfe)
}


#' Apply Custom QC Metric to Columns in a Spatial Transcriptomics Experiment
#'
#' This function adds a custom QC metric to columns in a spatial
#' transcriptomics experiment.
#'
#' @param sfe A SpatialFeatureExperiment object containing spatial
#' transcriptomics experiment data.
#' @param qcMetric A vector of TRUE and FALSE for columns (locations) that are
#' flagged for discarding.
#' @param qc_name The name of the custom QC metric to be added to the colData.
#'
#' @return An updated SpatialFeatureExperiment object with the custom QC metric
#' added to the colData.
#'
#' @details The function takes a vector `qcMetric` containing TRUE and FALSE
#' values, representing columns (spots) that should be flagged for discarding.
#' It adds a new column with the specified name `qc_name` in the colData of the
#' spatial transcriptomics experiment, marking spots based on the custom QC
#' metric.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords custom QC metric, SpatialFeatureExperiment, spatial features
#'
#' @rdname dot-int_custom_2
#'
.int_custom_2 <- function(sfe, qcMetric, qc_name) {
  colData(sfe)[[qc_name]] <- qcMetric

  return(sfe)
}
