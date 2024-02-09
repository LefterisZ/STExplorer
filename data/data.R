#' 10X Visium Steatotic Human Liver - processed
#'
#' A SpatialFeatureExperiment object containing a 10X Visium liver dataset from
#' the Liver Cell Atlas consortium. The sample used is the JBO019. The dataset
#' has passed QC for spot and feature selection.
#'
#' @format ## `sfe`
#' A SpatialFeatureExperiment S4 class object with 8535 elements. Specifically,
#' 8535 features and 1161 spots. Includes rowData with feature metadata and
#' statistics. Includes colData with spot metadata and statistics. Includes
#' colGeometry data with sf geometry objects. Includes raw and log2-normalised
#' UMI count data.
#'
#' The QC filters used:
#' \itemize{
#'    \item Spot-level QC
#'    \itemize{
#'       \item qc_lib_size: removed spots with library size (total UMI counts
#'       in the spot) of < 1000 or > 45000.
#'       \item qc_detected: removed spots with number of detected genes < 500
#'       or > 6000.
#'       \item qc_mito: removed spots with percentage of mitochondrial
#'       expression > 22%.
#'    }
#'    \item Gene-level QC
#'    \itemize{
#'       \item is_zero: removed any gene with no expression in the dataset.
#'       \item is_mito: removed all mitochondrial genes.
#'       \item is_logLow: remove genes with a sample mean of log2 UMI counts
#'       <1. Sample mean is the mean of counts over the locations a gene is
#'       present.
#'    }
#' }
#'
#' @usage data(sfe)
#'
#' @source <https://www.livercellatlas.org/download.php>
"sfe"

#' 10X Visium Steatotic Human Liver - unprocessed
#'
#' A SpatialFeatureExperiment object containing a 10X Visium liver dataset from
#' the Liver Cell Atlas consortium. The sample used is the JBO019. The dataset
#' has NOT passed any QC and any preprocessing.
#'
#' @format ## `sfe`
#' A SpatialFeatureExperiment S4 class object with 32738 elements. Specifically,
#' 8535 features and 1185 spots.
#'
#' @usage data(sfe_raw)
#'
#' @source <https://www.livercellatlas.org/download.php>
"sfe_raw"

#' Geographically Weighted Principal Components Analysis - GWPCA
#'
#' Results from running the [gwpcaSTE()] function on the processed and QC'ed
#' 10X Visium Steatotic Human Liver sample JBO019 from the Liver Cell Atlas
#' <https://www.livercellatlas.org/download.php>.
#'
#' @format ## `gwpca`
#' A list object of class `gwpca`. Includes all the data generated from running
#' the [gwpcaSTE()] function. Specifically includes:
#' \itemize{
#'    \item pca: global PCA results.
#'    \item loadings: loading scores of each genes in each Principal Component
#'    (PC) in each location.
#'    \item SDF a `SpatialPointsDataFrame` with Percentage of Variation
#'    explained by each PC in each location.
#'    \item gwpca.scores: PCA scores for each component in each location.
#'    \item var:
#'    \item local.PV: Percent of Variation in each location.
#'    \item GW.arguments: The arguments used to run [gwpcaSTE()].
#'    \item CV: Cross Validation results.
#'    \item geometry: Hexagonal geometries for each location. Used for plotting.
#' }
#'
#' The arguments used for this dataset are the below:
#' \itemize{
#'    \item vars: 526 genes vector as generated using the [scran::getTopHVGs()]
#'    function after running the [scran::modelGeneVar()] function on the
#'    preprocessed SFE object from the JBO019 sample. For the SFE object and
#'    the parameters used in the QC see also the `?STExplorerDev::sfe` help
#'    page. The parameters used in [scran::getTopHVGs()] are the below:
#'    \itemize{
#'       \item var.field: "bio"
#'       \item prop: 0.5
#'       \item fdr.threshold: 0.05
#'    }
#'    \item k: 20
#'    \item kernel: "gaussian"
#'    \item adaptive: FALSE
#'    \item p: 2
#'    \item cv: TRUE
#'    \item scores: FALSE
#'    \item robust: FALSE
#'    \item sfe: sfe object preprocessed and QC'ed.
#'    \item assay: "logcounts"
#'    \item bw: 882.09951902194
#' }
#'
#' @usage data(gwpca)
#'
"gwpca"

#' 10X Visium Steatotic Human Liver - annotation
#'
#' Annotation data from the JBO019 Human Liver dataset.
#'
#' @format ## `matrix`
#' A 3-column matrix with annotation data as annotated by the Liver Cell Atlas
#' project. Includes the below columns:
#' \itemize{
#'    \item Barcode: Spot barcode.
#'    \item sample_id: Sample ID.
#'    \item annotation: The manual annotation of each spot.
#'  }
#' Each column is `character`.
#'
#' @usage data(gTruth)
#'
#' @source <https://www.livercellatlas.org/download.php>
"gTruth"

#' Set of liver marker genes
#'
#' Marker genes for liver cell types and subtypes. Used for the heatmaps and
#' subclustering of FGWC
#'
#' @format ## `data.frame`
#' A 4-column data frame with marker genes and the cell types. Includes the
#' below columns:
#' \itemize{
#'    \item gene.name: the gene names
#'    \item ensg.ID: the ENSGene IDs
#'    \item Type: the cell type
#'    \item Subtype: the cell subtype
#' }
#' Each column is `character`.
#'
#' @usage data(markers)
#'
"markers"
