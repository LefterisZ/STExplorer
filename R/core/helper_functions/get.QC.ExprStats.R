#' @name get.QC.ExprStats
#' 
#' @description 
#' A function to calculate some more expression stats.
#'              
#' @param obj a SpatialFeaturesExperiment objects.
#' 
#' @param assay the name of the assay to use.
#' 
#' @param .sample_id the sample id to run the calculations for
#' 
#' @export

get.QC.ExprStats <- function(obj, assay, .sample_id) {
  
  ## Fetch the data per sample and perform the calculations
  data <- assay(obj, assay)[,colData(obj)$sample_id == .sample_id]
  
  ## Get the number of locations a gene is present
  s_nLocations <- Matrix::rowSums(data != 0)
  p_nLocations <- sum(colData(obj)$sample_id == .sample_id)
  
  ## Get the mean of UMI counts over the locations a gene is present
  s_mean <- rowData(obj)[[paste0(.sample_id, ".total")]] / s_nLocations
  ## Get the mean of UMI counts over the whole sample
  p_mean <- rowData(obj)[[paste0(.sample_id, ".total")]] / p_nLocations
  
  ## Get the max of UMI counts of the locations a gene is present
  max <- sparseMatrixStats::rowMaxs(data)
  
  ## Get the median of UMI counts over the locations a gene is present
  data <- Matrix::as.matrix(data)
  s_median <- apply(data, 1, 
                    function(dt) { 
                      dt <- dt[dt != 0]
                      median <- median(dt)
                    })
  ## Get the median of UMI counts over the whole sample
  p_median <- apply(data, 1, median)
  
  ## Get the min of UMI counts of the locations a gene is present
  s_min <- apply(data, 1, 
                 function(dt) {
                   dt <- dt[dt != 0]
                   min <- min(dt)
                 })
  
  ## Get the Stand. Dev. of UMI counts of the locations a gene is present
  s_SD <- apply(data, 1, 
                function(dt) {
                  dt <- dt[dt != 0]
                  sd <- sd(dt)
                })
  
  ## Get the Stand. Dev. of UMI counts of the whole sample
  p_SD <- apply(data, 1, sd)
  
  ## Add them to the rowData
  rowData(obj)[paste0(.sample_id, ".nLocations")] <- s_nLocations
  rowData(obj)[paste0(.sample_id, ".s_min")] <- s_min
  rowData(obj)[paste0(.sample_id, ".max")] <- max
  rowData(obj)[paste0(.sample_id, ".s_mean")] <- s_mean
  rowData(obj)[paste0(.sample_id, ".s_median")] <- s_median
  rowData(obj)[paste0(.sample_id, ".s_SD")] <- s_SD
  rowData(obj)[paste0(.sample_id, ".p_mean")] <- p_mean
  rowData(obj)[paste0(.sample_id, ".p_median")] <- p_median
  rowData(obj)[paste0(.sample_id, ".p_SD")] <- p_SD
  
  return(obj)
}
