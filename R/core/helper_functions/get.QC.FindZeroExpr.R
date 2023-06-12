#' @name get.QC.FindZeroExpr
#' 
#' @description 
#' A function to find any genes with no expression in all spots by calculating
#' the total expression of a gene on all spots.
#' 
#' @param obj a SpatialFeaturesExperiment object.
#' 
#' @param assay the name of the assay to use.
#' 
#' @param .sample_id the sample id to run the calculations for.
#' 
#' @export

get.QC.FindZeroExpr <- function(obj, assay, .sample_id){
  
  ## Fetch the data per sample and perform the calculations
  data <- assay(obj, assay)[,colData(obj)$sample_id == .sample_id]
  
  ## Find the total in-sample number of UMIs 
  total <- Matrix::rowSums(data)
  
  rowData(obj)[paste0(.sample_id, ".total")] <- total
  
  return(obj)
}
