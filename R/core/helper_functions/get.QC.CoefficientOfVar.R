#' Calculate gene Coefficient of Variation
#'
#' @name get.QC.CoefficientOfVar
#'
#' @description
#' IMPORTANT: This function is not exported to be used alone.
#' A function to calculate the coefficient of variation (CV) for each gene.
#' CV for each gene is calculated as the ratio of the standard deviation to the
#' mean and is expressed as a percentage. A lower CV would therefore mean higher
#' stability.
#'
#' @param sfe A SpatialFeaturesExperiment object.
#' @param assay The name of the assay to use.
#' @param .sample_id The sample ID to calculate the CV for.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @details This function calculates the CV for each gene in the given
#' `SpatialFeaturesExperiment` object and assigns the CV values to the rowData
#' of the object. The CV is calculated separately for the sample (`s_CV`) and
#' population (`p_CV`) data based on the provided `assay` and `.sample_id`.
#'
#' @return
#' The modified SpatialFeaturesExperiment object with added Coefficient of
#' Variance statistic calculated.
#'
#' @importFrom SpatialFeatureExperiment rowData colData
#' @importFrom SummarizedExperiment assay
#'
#' @seealso \code{\link{SpatialFeatureExperiment}}
#'
get.QC.CoefficientOfVar <- function(sfe, assay, .sample_id) {

  ## Fetch the data per sample and perform the calculations
  data <- assay(sfe, assay)[,colData(sfe)$sample_id == .sample_id]

  ## Calculate CV only for locations the gene is present
  sd = rowData(sfe)[[paste0(.sample_id, ".s_SD")]]
  mean = rowData(sfe)[[paste0(.sample_id, ".s_mean")]]
  s_CV = (sd / mean) * 100

  ## Calculate CV only for all locations in the sample
  sd = rowData(sfe)[[paste0(.sample_id, ".p_SD")]]
  mean = rowData(sfe)[[paste0(.sample_id, ".p_mean")]]
  p_CV = (sd / mean) * 100

  rowData(sfe)[paste0(.sample_id, ".s_CV")] <- s_CV
  rowData(sfe)[paste0(.sample_id, ".p_CV")] <- p_CV

  return(sfe)
}
