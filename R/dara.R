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
