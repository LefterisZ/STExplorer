#' Add location-related QC metrics
#'
#' @name addPerLocQC
#'
#' @description
#' A function to add a series of location (spot)-related QC metrics..
#'
#' @param sfe The SpatialFeatureExperiment object.
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
#' @param MARGIN 1 = Features, 2 = Locations. Indicates for which aspect the
#' sparsity will be calculated for; features (genes) or locations.
#'
#' @param ... further arguments passed to \code{addPerCellQC}, to pass to
#' \code{perCellQCMetrics.}
#'
#'
#' @importFrom scater addPerCellQC
#'
#' @examples
#' # Load data
#' data(sfe_raw)
#' data(gTruth)
#'
#' # Calculate location QC stats
#'
#' sfe <- addPerLocQC(sfe_raw, gTruth = gTruth, assay = "counts", MARGIN = 2)
#'
#'
#' @seealso \code{\link{addPerCellQC}}, \code{\link{add.barcodes}},
#' \code{\link{add.gTruth}}, \code{\link{add.index}},
#' \code{\link{get.QC.Sparsity}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @export
addPerLocQC <- function(sfe,
                        gTruth = NULL,
                        assay = "counts",
                        in_tissue = TRUE,
                        MARGIN,
                        ...) {

  ## Add Barcodes as column
  sfe <- add.barcodes(sfe)

  ## Add custom annotations (if are available)
  if (!is.null(gTruth)) {
    sfe <- add.gTruth(sfe, gtruth = gTruth)
  }

  ## Add index column
  sfe <- add.index(sfe)

  ## Add locational sparsity
  sfe <- get.QC.Sparsity(sfe, assay = assay, MARGIN = MARGIN)

  ## Add other locational QC metrics from scatter package
  sfe <- addPerCellQC(sfe, ...)

  return(sfe)

}
