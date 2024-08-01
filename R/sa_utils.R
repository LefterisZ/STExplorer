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

  if (test == "permutation") {
    colPVal <- paste0(statistic, "Pval_perm")
  } else if (test == "permutation") {
    colPVal <- paste0(statistic, "Pval_test")
  }

  if (statistic == "moran") {
    colStat <- "moranI"
  } else if (statistic == "geary") {
    colStat <- "gearyC"
  } else if (statistic == "getis") {
    colStat <- "getisG"
  }

  ## Export data from rowData
  data <- rowData(sfe)[c(colStat, colPVal)] %>%
    as.data.frame()

  ## Find the genes pass the pVal threshold
  if (!is.null(pVal)) {
    data <- mutate(data, filter_PVal = .data[[colPVal]] < pVal)
  }

  ## Find the genes that pass the statistic thershold
  if (!is.null(stat_thresh)) {
    data <- mutate(data, filter_Stat = abs(.data[[colStat]]) > stat_thresh)
  }

  ## Select the filter columns
  selected_cols <- grepl("filter_", colnames(data))

  ## Create a vector where each element corresponds to a row in the data frame
  keep <- apply(data[selected_cols], 1, function(row) all(row))
  names(keep) <- NULL

  ## Make sure the result has only logical (TRUE/FALSE) values
  keep <- as.logical(keep)

  ## Add the vector as a column in data
  data$keep <- keep

  ## Return a named vector
  out <- rownames(data)[data$keep]
  names(out) <- out

  return(out)
}
