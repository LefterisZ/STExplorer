#' A set of functions to calculate the sparsity of the dataset and get the 
#' sparsity of each location.
#'              

#' @name .get.Global.Sparsity
#' @description
#' A function to fetch the dataset sparsity attributes from an assay.
#' 
.get.Global.Sparsity <- function(obj) {
  
  tmp <- metadata(obj)$Sparsity
  
  tmp <- dplyr::bind_rows(tmp)
  
  colnames(tmp) <- c("Assay", "Sample ID", "Total", "Zeros", "Sparsity %")
  
  return(tmp)
}

#' @name .int_sparsity_1
#' 
#' @description 
#' A function to calculate the sparsity of the dataset and get the 
#' sparsity of each gene.
#'              
#' @param obj a SpatialFeaturesExperiment objects.
#' 
#' @param assay the name of the assay to use.
#' 
#' @param sampleNo used to call the sparsity also over aall datasets when 
#' multiple samples exist.Used only when MARGIN = 1 at \code{addPerGeneQC} 
#' function.
#' 
#' @param .sample_id the sample id to run the calculations for. Used only when
#' MARGIN = 1 at \code{addPerGeneQC} function.
#' 

.int_sparsity_1 <- function(obj, 
                            assay,
                            sampleNo,
                            .sample_id) {
  ## Fetch data per sample
  data <- assay(obj, assay)[,colData(obj)$sample_id == .sample_id]
  
  ## Calculate sparsity
  ## per feature over all samples
  if (sampleNo > 1) {
    tot_zeros <- rowSums(assay(obj, assay) == 0)
    tot_sparsity <- tot_zeros/ncol(assay(obj, assay))
  }
  ## per feature per sample
  zeros <- rowSums(data == 0)
  sparsity <- zeros/ncol(data)
  
  ## Add the result to the output matrix
  ## per feature
  if (sampleNo > 1) {
    rowData(obj)$sparsity_tot <- tot_sparsity
  }
  ## per feature per sample
  rowData(obj)[paste0(.sample_id, ".sparsity")] <- sparsity
  
  ## Calculate dataset's sparsity of each location (per sample)
  zeros <- sum(data == 0)
  total <- dim(data)[1] * dim(data)[2]
  sparsity <- round(zeros/total*100, 2)
  metaSpar <- data.frame(assay, .sample_id, zeros, total, sparsity)
  
  ## Add dataset-wide, per locations, sparsity stats
  metadata(obj)$Sparsity[[paste0(assay, .sample_id)]] <- metaSpar
  
  return(obj)
}

#' @name .int_sparsity_2
#' 
#' @description 
#' A function to calculate the sparsity of the dataset and get the 
#' sparsity of each location.
#'              
#' @param obj a SpatialFeaturesExperiment objects.
#' 
#' @param assay the name of the assay to use.
#' 

.int_sparsity_2 <- function(obj, assay) {
  ## Fetch data per sample
  data <- assay(obj, assay)
  
  ## Calculate sparsity
  ## per location
  zeros <- colSums(data == 0)
  sparsity <- zeros/nrow(data)
  
  
  ## Add the result to the output matrix
  ## per location
  colData(obj)$sparsity <- sparsity
  
  return(obj)
}

#' @name get.QC.Sparsity
#' 
#' @description 
#' A function to calculate the sparsity of the dataset and get the 
#' sparsity of each location or each gene.
#'              
#' @param obj a SpatialFeaturesExperiment objects.
#' 
#' @param assay the name of the assay to use.
#' 
#' @param MARGIN 1 = Features, 2 = Locations. Indicates for which aspect the 
#' sparsity will be calculated; for features (genes) or locations.
#' 
#' @param sampleNo used to call the sparsity also over aall datasets when 
#' multiple samples exist.Used only when the function is called by 
#' \code{addPerGeneQC}.
#' 
#' @param .sample_id the sample id to run the calculations for. Used only when
#' the function is called by \code{addPerGeneQC}
#' 
#' @export

get.QC.Sparsity <- function(obj, 
                            assay, 
                            MARGIN,
                            sampleNo = NULL,
                            .sample_id = NULL) {
  
  ## Fetch data per sample
  if (MARGIN == 1) {
    obj <- .int_sparsity_1(obj = obj, assay = assay, 
                           sampleNo = sampleNo, .sample_id = .sample_id)
  } else if (MARGIN == 2) {
    obj <- .int_sparsity_2(obj = obj, assay = assay)
  }
  
  return(obj)
}

