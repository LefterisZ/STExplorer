#' @name get.QC.SA
#' 
#' @description
#' A function that calculates global Spatial Autocorrelation (SA) statistics 
#' for each gene and adds them in the \code{rowData} of the SFE 
#' (\code{SpatialFeatureExperiment}) object.
#' 
#' @param obj The SFE object.
#' @param type The different types that you can use to generate a neighbours 
#' graph. There are three categories; (a) Contiguity-based - "poly2nb", 
#' (b) Graph-based - "tri2nb", "soi.graph", "gabrielneigh", "relativeneigh" and
#' (c) Distance-based - "knearneigh", "dnearneigh". For more information about 
#' the individual functions visit the \code{spdep}'s documentation and vignette.
#' @param sample_id character string specifying unique sample identifier for 
#' the sample you need to generate a neighbour graph. Leave NULL if you want 
#' the same type of graph for all samples in the SFE object
#' @param ... arguments that are passed down to the spdep functions called by 
#' the \code{type} argument. To see what else is needed for the functions to 
#' operate correctly visit the \code{spdep}'s documentation and vignette..
#' 
#' @export

get.QC.SA <- function(obj, 
                      sample_id = NULL,
                      zero.policy = TRUE,
                      ...) {
  
  moran.test(x = counts(obj)[4,],
             listw = colGraph(obj, sample_id),
             randomisation = TRUE,
             zero.policy = zero.policy, 
             alternative = "greater", 
             rank = FALSE, 
             na.action = na.fail, 
             spChk = NULL,
             adjust.n = TRUE, 
             drop.EI2 = FALSE)
}
