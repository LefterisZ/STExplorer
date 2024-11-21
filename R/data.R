#' Fibrotic Lung data
#'
#' The fibrotic sample (V19S23-092-D1) from Franzen et. al., 2024, Nat. Genet.
#' publication <https://doi.org/10.1038/s41588-024-01819-2>.
#'
#' @format ## `msfe`
#' A MetaSpatialFeatureExperiment class object containing with 2 slots, one for
#' the SpatialFeatureExperiment (SFE) objects and one with the sample names for
#' each SFE object, containing one SFE object:
#' \describe{
#'   \item{sfe_data}{List of 1: contains the SFE object named "Fibrotic"}
#'   \item{sample_ids}{List of 1: contains the SFE object's name}
#' }
#' No preprocessing has been performed.
#'
#' @source <https://doi.org/10.1038/s41588-024-01819-2>
"msfe"


#' Healthy and Fibrotic Lung data
#'
#' The lung dataset comes from the Franzen et. al., 2024, Nat. Genet.
#' publication <https://doi.org/10.1038/s41588-024-01819-2>. The Healthy and
#' Fibrotic samples here are samples V10T31-015-A1 and V19S23-092-D1
#' respectively from the publication.
#'
#' @format ## `msfe_2`
#' A MetaSpatialFeatureExperiment class object containing with 2 slots, one for
#' the SpatialFeatureExperiment (SFE) objects and one with the sample names for
#' each SFE object, containing two SFE objects:
#' \describe{
#'   \item{sfe_data}{List of 1: contains the SFE object named "Fibrotic"}
#'   \item{sample_ids}{List of 1: contains the SFE object's name}
#' }
#' No preprocessing has been performed.
#'
#' @source <https://doi.org/10.1038/s41588-024-01819-2>
"msfe_2"


#' Marker genes for IPF lung cell types
#'
#' The information comes from the Habermann et. al., 2020, Sci. Adv.
#' publication <https://doi.org/10.1126/sciadv.aba1972>. The information can be
#' found in Supplementary Table 3 (Table S3).
#'
#' @format ## `markers`
#' A dataframe with cell-type markers containing the below columns
#' \describe{
#'   \item{gene.name}{List of 1: contains the SFE object named "Fibrotic"}
#'   \item{ensg.ID}{List of 1: contains the SFE object's name}
#'   \item{Type}{High level classes: Immune, Epithelieal, Mesenchymal, and
#'   Endothelial}
#'   \item{Subtype}{Low level classification of cells (i.e., HAS1-high
#'   Fibroblasts)}
#'   \item{TypeClass}{Mid level classification of cells (i.e., HAS1-high
#'   Fibroblasts and Fibroblasts both belong to Mid-level class Fibroblasts)}
#'   \item{MarkerInfo}{markers are divided in two categories: positive and
#'   negative}
#'   \item{ExpressionLevelInfo}{Some markers are marked as high or low
#'   expression levels}
#'   \item{description}{The gene's full name}
#' }
#' No preprocessing has been performed.
#'
#' @source <https://doi.org/10.1126/sciadv.aba1972>
"markers"


#' Subset of Fibrotic Lung data
#'
#' The lung dataset comes from the Franzen et. al., 2024, Nat. Genet.
#' publication <https://doi.org/10.1038/s41588-024-01819-2>. The Healthy and
#' Fibrotic samples here is the sample V19S23-092-D1 from the publication.
#'
#' This object contains a subset of the Fibrotic sample. Specifically the top
#' left corner as defined by the below selection:
#'
#' spatialCoords(sfe)[, "pxl_row_in_fullres"] > 8000 & spatialCoords(sfe)[, "pxl_col_in_fullres"] < 7000
#'
#' It is done so that the vignettes render faster and to be used for the
#' function examples.
#'
#' @format ## `sfe`
#' A SpatialFeatureExperiment class object containing the Fibrotic sample subset.
#'
#' @source <https://doi.org/10.1038/s41588-024-01819-2>
"sfe"


#' NMF for the subset of Fibrotic Lung data
#'
#' The lung dataset comes from the Franzen et. al., 2024, Nat. Genet.
#' publication <https://doi.org/10.1038/s41588-024-01819-2>. The Healthy and
#' Fibrotic samples here is the sample V19S23-092-D1 from the publication.
#'
#' This object contains the results of running NMF dimensionality reduction to
#' be used for FGWC.
#'
#' It is generated using the below command:
#' sfe_nmf <- fgwc_nmf(sfe, sample_id = "Fibrotic",
#' top_hvgs = top_hvgs[["Fibrotic"]], ncomponents = 7)
#'
#' To find the top_hvgs we used the whole Fibrotic sample data and the default
#' parameters in modeling gene variance and selecting topHVGs using
#' STExplorer's functions.
#'
#' @format ## `sfe_nmf`
#' A dataframe containing the NMF result for the Fibrotic sample subset.
#'
#' @source <https://doi.org/10.1038/s41588-024-01819-2>
"sfe_nmf"
