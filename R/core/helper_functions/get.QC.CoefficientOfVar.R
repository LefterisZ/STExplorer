#' @name get.QC.CoefficientOfVar
#' 
#' @description 
#' A function to calculate the coefficient of variation (CV) for each 
#' gene. CV for each gene is calculated as the ratio of the standard deviation 
#' to the mean and it is expressed as a percentage. A lower CV would therefore 
#' mean higher stability.
#'              
#' @param obj a SpatialFeaturesExperiment objects.
#' 
#' @param assay the name of the assay to use.
#' 
#' @export

get.QC.CoefficientOfVar <- function(obj, assay, .sample_id) {
  
  ## Fetch the data per sample and perform the calculations
  data <- assay(obj, assay)[,colData(obj)$sample_id == .sample_id]
  
  ## Calculate CV only for locations the gene is present
  sd = rowData(obj)[[paste0(.sample_id, ".s_SD")]]
  mean = rowData(obj)[[paste0(.sample_id, ".s_mean")]]
  s_CV = (sd / mean) * 100
  
  ## Calculate CV only for all locations in the sample
  sd = rowData(obj)[[paste0(.sample_id, ".p_SD")]]
  mean = rowData(obj)[[paste0(.sample_id, ".p_mean")]]
  p_CV = (sd / mean) * 100
  
  rowData(obj)[paste0(.sample_id, ".s_CV")] <- s_CV
  rowData(obj)[paste0(.sample_id, ".p_CV")] <- p_CV
  
  return(obj)
}
