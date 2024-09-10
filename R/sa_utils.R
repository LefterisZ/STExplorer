#' Get Spatial Autocorrelation Global Genes
#'
#' This function calculates spatial autocorrelation statistics and returns a
#' vector of genes that pass specified thresholds.
#'
#' @param m_sfe An object of class 'SpatialFeatureExperiment' or
#'              'MultiSpatialFeatureExperiment' containing spatial data.
#' @param sample_id A character vector specifying the sample ID to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param statistic A character string specifying the type of spatial
#'                  autocorrelation results ("moran", "geary", "getis") to look
#'                  for in the m_sfe object for the specified sample.
#' @param test A character string specifying the statistical test method
#'              ("z-score" or "permutation") to look for in the m_sfe object
#'              for the specified sample.
#' @param pVal A numeric value representing the significance threshold for
#'             plotting. Ignored for Geary's C, where permutation is suggested.
#'             See details for more information.
#' @param stat_thresh Numeric value specifying the threshold for the spatial
#'                    autocorrelation statistic. Genes with an absolute value
#'                    of the statistic greater than this threshold will be
#'                    included. Default is 0.5.
#'
#' @return A named vector of gene names that pass the specified thresholds.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial, autocorrelation
#'
#' @rdname getSAGlobalGenes
#'
#' @importFrom dplyr mutate
#'
#' @examples
#' getSAGlobalGenes()
#'
#' @export
getSAGlobalGenes <- function(m_sfe,
                             sample_id = NULL,
                             statistic = c("moran", "geary", "getis"),
                             test = c("z-score", "permutation"),
                             pVal = 0.05,
                             stat_thresh = 0.5) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get the column names based on the statistic and test
  columns <- .int_getColumnNames(statistic, test)
  colStat <- columns$colStat
  colPVal <- columns$colPVal

  ## Check that the necessary columns exist in the rowData colnames
  #required_columns <- c(colStat, colPVal)
  missing_columns <- setdiff(unlist(columns), colnames(rowData(sfe)))
  if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  ## Extract and prepare data
  data <- as.data.frame(rowData(sfe)[, c(colStat, colPVal)])

  ## Filter genes based on the provided thresholds
  filtered_data <- .int_filterGenes(data, colStat, colPVal, pVal, stat_thresh)

  ## Return the row names of the filtered data
  return(rownames(filtered_data))
}


# ---------------------------------------------------------------------------- #
#  ###### INTERNAL FUNCTIONS ASSOCIATED WITH SA CALCULATIONS (C, G, I) ######
# ---------------------------------------------------------------------------- #
#' Internal Function: Construct and export column names for Global SA results
#'
#' This function aims to construct column names based on the statistic and test.
#'
#' @inheritParams getSAGlobalGenes
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getColumnNames
#'
.int_getColumnNames <- function(statistic, test) {
  colPVal <- NULL
  colStat <- NULL

  if (test %in% c("permutation", "z-score")) {
    colPVal <- paste(statistic,
                     ifelse(test == "permutation", "PvalPerm", "PvalTest"),
                     sep = "_")
  } else {
    stop("Unsupported test type: ", test)
  }

  if (statistic == "moran") {
    colStat <- "moranI_stat"
    colPVal <- gsub("_", "I_", colPVal)
  } else if (statistic == "geary") {
    colStat <- "gearyC_stat"
    colPVal <- gsub("_", "C_", colPVal)
  } else if (statistic == "getis") {
    colStat <- "getisG_stat"
    colPVal <- gsub("_", "G_", colPVal)
  } else {
    stop("Unsupported statistic: ", statistic)
  }

  return(list(colStat = colStat, colPVal = colPVal))
}

#' Internal Function: filter genes based on Global SA statistic
#'
#' This function uses p-value and/or SA statistic values as thresholds to
#' filter and select genes. Returns a data frame of selected genes from the
#' SFE object's rowData.
#'
#' @param data the SFE rowData.
#' @param colStat the SA statistic's values column name.
#' @param colPVal the SA statistic's p-values column name.
#' @param pVal the SA statistic's p-value threshold. MUST be a numeric.
#' @param stat_thresh the SA statistic threshold. MUST be a numeric.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_filterGenes
.int_filterGenes <- function(data, colStat, colPVal, pVal, stat_thresh) {
  ## Validate pVal
  if (!is.null(pVal) && (!is.numeric(pVal) || pVal < 0 || pVal > 1)) {
    stop("pVal should be a numeric value between 0 and 1.")
  }

  ## Validate stat_thresh
  if (!is.null(stat_thresh) && (!is.numeric(stat_thresh))) {
    stop("stat_thresh should be a numeric value.")
  }

  ## Apply the p-value threshold
  if (!is.null(pVal)) {
    if (any(is.na(data[[colPVal]]))) {
      warning("Missing values found in p-value column: ", colPVal)
    }
    data <- mutate(data, filter_PVal = .data[[colPVal]] < pVal)
  }

  ## Apply the statistic threshold
  if (!is.null(stat_thresh)) {
    if (any(is.na(data[[colStat]]))) {
      warning("Missing values found in statistic column: ", colStat)
    }
    data <- mutate(data, filter_Stat = abs(.data[[colStat]]) > stat_thresh)
  }

  ## Combine the filter columns
  selected_cols <- grep("filter_", names(data), value = TRUE)

  ## Ensure that we have filter columns to process
  if (length(selected_cols) == 0) {
    stop("No filtering has been applied. Check threshold parameters.")
  }

  ## Create a logical vector to keep rows that pass all filters
  keep <- rowSums(data[, selected_cols]) == length(selected_cols)
  names(keep) <- NULL

  ## Ensure the result has only logical (TRUE/FALSE) values
  keep[is.na(keep)] <- FALSE

  return(data[keep, ])
}


# ---------------------------------------------------------------------------- #
#  ###### INTERNAL FUNCTIONS ASSOCIATED WITH SA CALCULATIONS (C, G, I) ######
# ---------------------------------------------------------------------------- #
#' Internal Function: check genes input
#'
#' This internal function checks for the genes input in the main function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @importFrom SummarizedExperiment assay
#'
#' @rdname dot-int_SAgeneCheck
#'
#' @aliases dot-int_SAgeneCheck
.int_SAgeneCheck <- function(sfe, genes) {
  ## Check genes to use
  if (is.character(genes)) {
    # keep only genes that are present in the SFE object.
    # Otherwise throws an error downstream.
    gs <- genes %in% rownames(sfe)
    if (length(genes) > sum(gs)) {
      message("These genes are not present in your dataset: \n\t",
              "ENSGIDs: ", paste(genes[!gs], collapse = " "), "\n\t",
              "Gene Names: ", paste(names(genes[!gs]), collapse = " "))
    }
    genes <- genes[gs]
    names(genes) <- genes
  } else if (isTRUE(genes)) {
    genes <- rownames(rowData(sfe))
    names(genes) <- genes
  } else {
    stop("Invalid `genes` argument input")
  }

  return(genes)
}
