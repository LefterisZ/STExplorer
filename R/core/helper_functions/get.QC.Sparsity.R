#' Calculate sparsity
#'
#' @name get.QC.Sparsity
#'
#' @description
#' A function to calculate the sparsity of the dataset and get the
#' sparsity of each location or each gene.
#'
#' @param sfe A SpatialFeaturesExperiment object.
#'
#' @param assay The name of the assay to use.
#'
#' @param MARGIN Specifies the aspect for which the sparsity will be calculated.
#' Use 1 for features (genes) or 2 for locations.
#'
#' @param sampleNo The sample number used to call the sparsity over all datasets
#' when multiple samples exist. This parameter is used only when the function
#' is called by addPerGeneQC.
#'
#' @param .sample_id The sample ID to run the calculations for. This parameter
#' is used only when the function is called by addPerGeneQC.
#'
#' @details The sparsity is calculated as the proportion of zeros in the dataset.
#' It provides an indication of the density of expression values. A higher
#' sparsity indicates a higher number of zero values, indicating low expression
#' levels or absence of expression.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso addPerGeneQC
#'
get.QC.Sparsity <- function(sfe,
                            assay,
                            MARGIN,
                            sampleNo = NULL,
                            .sample_id = NULL) {

  ## Fetch data per sample
  if (MARGIN == 1) {
    sfe <- .int_sparsity_1(sfe = sfe, assay = assay,
                           sampleNo = sampleNo, .sample_id = .sample_id)
  } else if (MARGIN == 2) {
    sfe <- .int_sparsity_2(sfe = sfe, assay = assay)
  }

  return(sfe)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Fetch global sparsity stats
#'
#' @name .get.Global.Sparsity
#' @description
#' A function to fetch the dataset sparsity attributes from an assay.
#' @param sfe A SpatialFeatureExperiment object
#'
#' @importFrom S4Vectors metadata
#' @importFrom dplyr bind_rows
#'
.get.Global.Sparsity <- function(sfe) {

  tmp <- metadata(sfe)$Sparsity

  tmp <- dplyr::bind_rows(tmp)

  colnames(tmp) <- c("Assay", "Sample ID", "Total", "Zeros", "Sparsity %")

  return(tmp)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Internal: calculate feature sparsity
#'
#' @name .int_sparsity_1
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to calculate the sparsity of the dataset and get the
#' sparsity of each gene.
#'
#' @param sfe a SpatialFeaturesExperiment objects.
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
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialFeatureExperiment colData rowData
#' @importFrom S4Vectors metadata metadata<-
#'
.int_sparsity_1 <- function(sfe,
                            assay,
                            sampleNo,
                            .sample_id) {
  ## Fetch data per sample
  data <- assay(sfe, assay)[,colData(sfe)$sample_id == .sample_id]

  ## Calculate sparsity
  ## per feature over all samples
  if (sampleNo > 1) {
    tot_zeros <- rowSums(assay(sfe, assay) == 0)
    tot_sparsity <- tot_zeros/ncol(assay(sfe, assay))
  }
  ## per feature per sample
  zeros <- rowSums(data == 0)
  sparsity <- zeros/ncol(data)

  ## Add the result to the output matrix
  ## per feature
  if (sampleNo > 1) {
    rowData(sfe)$sparsity_tot <- tot_sparsity
  }
  ## per feature per sample
  rowData(sfe)[paste0(.sample_id, ".sparsity")] <- sparsity

  ## Calculate dataset's sparsity of each location (per sample)
  zeros <- sum(data == 0)
  total <- dim(data)[1] * dim(data)[2]
  sparsity <- round(zeros/total*100, 2)
  metaSpar <- data.frame(assay, .sample_id, zeros, total, sparsity)

  ## Add dataset-wide, per locations, sparsity stats
  metadata(sfe)$Sparsity[[paste0(assay, .sample_id)]] <- metaSpar

  return(sfe)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Internal: calculate location sparsity
#'
#' @name .int_sparsity_2
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to calculate the sparsity of the dataset and get the
#' sparsity of each location.
#'
#' @param sfe a SpatialFeaturesExperiment objects.
#'
#' @param assay the name of the assay to use.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialFeatureExperiment colData
#'
.int_sparsity_2 <- function(sfe, assay) {
  ## Fetch data per sample
  data <- assay(sfe, assay)

  ## Calculate sparsity
  ## per location
  zeros <- colSums(data == 0)
  sparsity <- zeros/nrow(data)


  ## Add the result to the output matrix
  ## per location
  colData(sfe)$sparsity <- sparsity

  return(sfe)
}
