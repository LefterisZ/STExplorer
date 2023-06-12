#' @name addPerGeneQC
#' 
#' @description
#' A function to add a series of location (spot)-related QC metrics..
#' 
#' @param obj The SpatialFeaturesExperiment object.
#' 
#' @param gTruth A dataframe that contains the ground truth for your dataset. 
#' It needs to have at least 3 columns. One column named "Barcode" with the 
#' spot Barcodes (these need to match the colnames of the SFE object), one 
#' column named "sample_id" which includes the sample ID for each spot (this is
#' important when you import multiple samples) and one column with the ground 
#' truth annotation.
#' 
#' @param assay the name of the assay to use. Defaults to 'counts'.
#' 
#' @param version The ENSEMBL version of the annotation you want to use. It is 
#' advised to use the annotation version that was used to create the .bam files.
#' If you don't want to add annotation to the genes leave it as NULL.
#' 
#' @param ... further arguments passed to \code{addPerFeatureQCMetrics}, to pass
#'  to \code{perFeatureQCMetrics}.
#' 
#' @export

addPerGeneQC <- function(obj,
                        assay = "counts",
                        organism = "human",
                        version = NULL,
                        ...) {
  ## Add Biomart annotations
  if (!is.null(version)) {
    mart <- create_biomart(organism = organism, version = version)
    rowData(obj) <- annotate_data.frame(rowData(obj), biomart = mart)
  }
  
  ## Add gene QC metrics from scatter package
  obj <- addPerFeatureQC(obj, ...)
  
  ## Find the total inter-sample number of UMIs 
  total <- Matrix::rowSums(assay(obj, assay))
  rowData(obj)$total <- total
  
  ## Prepare for putatively multiple samples
  samples = unique(colData(obj)$sample_id)
  
  for (s in samples) {
    message(paste0("calculating stats for sample: ", s))
    ## Add locational sparsity
    obj <- get.QC.Sparsity(obj, assay = assay, MARGIN = 1, 
                           sampleNo = length(samples), .sample_id = s)
    
    ## Find genes with zero counts
    obj <- get.QC.FindZeroExpr(obj, assay = assay, .sample_id = s)
    
    ## Add other gene stats 
    obj <- get.QC.ExprStats(obj, assay = assay, .sample_id = s)

    ## Add Coefficient of Variance
    obj <- get.QC.CoefficientOfVar(obj, assay = assay, .sample_id = s)
  }
  
  return(obj)
}