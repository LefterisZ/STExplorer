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
