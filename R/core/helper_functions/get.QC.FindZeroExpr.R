#' Find non-expressed genes
#'
#' @name get.QC.FindZeroExpr
#'
#' @description
#' A function to find zero expression values and calculate total UMIs.
#'
#' @param sfe A SpatialFeatureExperiment object.
#'
#' @param assay The name of the assay to use.
#'
#' @param .sample_id The sample ID to run the calculations for.
#'
#' @details This function takes a SpatialFeaturesExperiment object and
#' calculates the total number of UMIs for each gene in the specified sample.
#' It extracts the data for the given sample from the specified assay, finds the
#' total UMIs for each gene, and stores the results in the rowData of the
#' SpatialFeaturesExperiment object.
#'
#' @importFrom Matrix rowSums
#' @importFrom SpatialFeatureExperiment rowData
#' @importFrom SummarizedExperiment assay
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso \code{\link{SpatialFeatureExperiment}}, \code{\link{rowData}},
#' \code{\link{assay}}
#'
#' @return The modified SpatialFeaturesExperiment object with the total UMIs
#' stored in the rowData.
#'

get.QC.FindZeroExpr <- function(sfe, assay, .sample_id){

  ## Fetch the data per sample and perform the calculations
  data <- assay(sfe, assay)[,colData(sfe)$sample_id == .sample_id]

  ## Find the total in-sample number of UMIs
  total <- Matrix::rowSums(data)

  rowData(sfe)[paste0(.sample_id, ".total")] <- total

  return(sfe)
}
