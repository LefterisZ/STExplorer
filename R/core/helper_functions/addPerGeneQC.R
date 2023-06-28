#' Add feature-related QC metrics
#'
#' @name addPerGeneQC
#'
#' @description
#' A function to add a series of feature (gene)-related QC metrics.
#'
#' @param sfe The SpatialFeaturesExperiment object.
#'
#' @param assay the name of the assay to use. Defaults to 'counts'.
#'
#' @param organism The organism for which the annotation is used. Defaults to
#' "human". Other names that can be used: mouse, zebrafish, rat, fruitfly, pig,
#' worm, yeast.
#'
#' @param version The ENSEMBL version of the annotation you want to use. It is
#' advised to use the annotation version that was used to create the .bam files.
#' If you don't want to add annotation to the genes leave it as NULL.
#'
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here
#' are 'www', 'uswest', 'useast', 'asia'. If no mirror is specified the primary
#' site at www.ensembl.org will be used. Mirrors are not available for the
#' Ensembl Genomes databases. The mirror and the version arguments cannot be
#' used together.
#'
#' @param ... further arguments passed to \code{addPerFeatureQC}, to pass
#'  to \code{perFeatureQCMetrics}.
#'
#' @importFrom scater addPerFeatureQC
#' @importFrom Matrix rowSums
#' @importFrom SpatialFeatureExperiment rowData colData
#' @importFrom SummarizedExperiment rowData<-
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' data(sfe_raw)
#' sfe <- addPerGeneQC(sfe_raw, assay = "counts", organism = "human",
#' version = NULL, mirror = NULL)
#'
#' @seealso \code{\link{addPerFeatureQC}}, \code{\link{get.QC.Sparsity}},
#' \code{\link{get.QC.FindZeroExpr}}, \code{\link{get.QC.ExprStats}},
#' \code{\link{get.QC.CoefficientOfVar}}
#'
#' @export
#'
addPerGeneQC <- function(sfe,
                        assay = "counts",
                        organism = "human",
                        version = NULL,
                        mirror = NULL,
                        ...) {
  ## Add Biomart annotations
  if (is.null(version)) {
      colnames(rowData(sfe)) <- "gene_name"
  } else {
    mart <- try(createBiomart(organism = organism,
                              version = version,
                              mirror = mirror))
    rowData(sfe) <- try(annotateDataFrame(rowData(sfe), biomart = mart))
  }

  ## Add gene QC metrics from scater package
  sfe <- addPerFeatureQC(sfe, ...)

  ## Find the total inter-sample number of UMIs
  total <- Matrix::rowSums(assay(sfe, assay))
  rowData(sfe)$total <- total

  ## Prepare for putatively multiple samples
  samples = unique(colData(sfe)$sample_id)

  for (s in samples) {
    message(paste0("calculating stats for sample: ", s))
    ## Add locational sparsity
    sfe <- get.QC.Sparsity(sfe, assay = assay, MARGIN = 1,
                           sampleNo = length(samples), .sample_id = s)

    ## Find genes with zero counts
    sfe <- get.QC.FindZeroExpr(sfe, assay = assay, .sample_id = s)

    ## Add other gene stats
    sfe <- get.QC.ExprStats(sfe, assay = assay, .sample_id = s)

    ## Add Coefficient of Variance
    sfe <- get.QC.CoefficientOfVar(sfe, assay = assay, .sample_id = s)
  }

  return(sfe)
}
