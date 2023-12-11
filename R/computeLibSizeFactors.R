#' Compute Library Size Factors in a SpatialFeatureExperiment
#'
#' This function computes library size factors in a
#' \code{SpatialFeatureExperiment} object. Size factors represent a scaling
#' factor for each sample, adjusting for differences in sequencing depth.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param type Character string, one of "inter" or "intra", indicating whether
#' to calculate size factors within each sample ("intra") or across all samples
#' ("inter"). The default is "intra".
#' @param ... Additional arguments passed to
#' \code{\link[scater]{computeLibraryFactors}}.
#'
#' @return A modified \code{SpatialFeatureExperiment} object with computed
#' library size factors.
#'
#' @details This function calculates library size factors for each sample
#' individually ("intra") or for all samples together ("inter"). The type of
#' calculation is determined by the \code{type} argument.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics library size factors
#'
#' @rdname computeLibSizeFactors
#'
#' @seealso \code{\link[scater]{computeLibraryFactors}}
#'
#' @importFrom scater librarySizeFactors
#'
#' @examples
#' \dontrun{# Compute library size factors within each sample
#' sfe # An SFE object
#' sfe <- computeLibSizeFactors(sfe, type = "intra")
#'
#' # Compute library size factors across all samples
#' sfe # An SFE object
#' sfe <- computeLibSizeFactors(sfe, type = "inter")
#' }
#'
#' @export
computeLibSizeFactors <- function(sfe, type = c("inter", "intra"), ...) {

  ## Check arguments
  stopifnot(is(sfe, "SpatialFeatureExperiment"))

  ## Check valid type argument
  if (missing(type)) {
    in_sample <- TRUE
  } else {
    type <- match.arg(type)
    if (type == "intra") {
      in_sample <- TRUE
    } else {
      in_sample <- FALSE
    }
  }

  if (!in_sample) {
    ## Calculate library size factors
    ## Take into account all samples together
    sfe <- scater::computeLibraryFactors(sfe, ...)
  } else if (in_sample) {
    ## Calculate library size factors
    ## Take into account one sample at a time

    ## Get sample IDs
    ids <- .int_getSmplIDs(sfe = sfe, sample_id = TRUE)

    ## Create a list of filtered SpatialFeatureExperiment objects
    sfe_list <- lapply(ids, .int_sizeFactCalc, sfe = sfe, ...)

    ## Merge the list of SpatialFeatureExperiment objects into one
    sfeOut <- Reduce(function(x, y) cbind(x, y), sfe_list)

    ## Clean up duplicate entries in the metadata slot
    sfeOut <- .int_cleanMetaData(sfeOut = sfeOut)
  }

  ## Have a look at the size factors
  print(.int_summarySizeFact(sfe = sfeOut, ids = ids))

  return(sfeOut)
}

#' Internal Function: Compute Library Size Factors for a Specific Sample
#'
#' Compute Library Size Factors for a Specific Sample
#'
#' @param id Character string, specifying the sample/image identifier.
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param ... Additional arguments passed to \code{\link[scater]{computeLibraryFactors}}.
#'
#' @return A modified \code{SpatialFeatureExperiment} object with computed library size factors.
#'
#' @details This function calculates library size factors for a specific sample identified by
#'           the \code{id} argument using the \code{\link[scater]{computeLibraryFactors}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics library size factors
#'
#' @rdname computeLibSizeFactors
#'
#' @seealso \code{\link[scater]{computeLibraryFactors}}
#'
.int_sizeFactCalc <- function(id, sfe, ...) {
  sfe_int <- sfe[, colData(sfe)$sample_id == id]
  sfe_int <- scater::computeLibraryFactors(sfe_int, ...)
  return(sfe_int)
}


#' Internal Function: Summary of Library Size Factors for a Specific Sample
#'
#' Summary of Library Size Factors for a Specific Sample
#'
#' @param id Character string, specifying the sample/image identifier.
#' @param sfe A \code{SpatialFeatureExperiment} object.
#'
#' @return A summary of library size factors for the specified sample.
#'
#' @details This function provides a summary of library size factors for a
#' specific sample identified by the \code{id} argument. It extracts the size
#' factors using the \code{\link[scater]{librarySizeFactors}} function and
#' summarises them.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics library size factors summary
#'
#' @rdname computeLibSizeFactors
#'
#' @seealso \code{\link[scater]{librarySizeFactors}}
#'
.int_summaryCalc <- function(id, sfe) {
  sfe_int <- sfe[, colData(sfe)$sample_id == id]
  sfe_int <- summary(librarySizeFactors(sfe_int))
  return(sfe_int)
}


#' Internal Function: Summary of Library Size Factors for Multiple Samples
#'
#' This function computes a summary of library size factors for multiple
#' samples in a spatial transcriptomics experiment.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param ids A character vector specifying sample/image identifiers.
#'
#' @return A data frame summarizing library size factors for multiple samples.
#'
#' @details The function internally uses the
#' \code{\link[scater]{librarySizeFactors}} function for calculating library
#' size factors and the \code{\link{.int_summaryCalc}} function to generate
#' summaries for each sample. The results are then combined into a data frame.
#'
#' This function efficiently computes a summary of library size factors for
#' multiple samples in a spatial transcriptomics experiment. It leverages the
#' capabilities of the internal function \code{.int_summaryCalc}, ensuring a
#' systematic and comprehensive overview of library size factors.
#'
#' The data is organized into a data frame with rows corresponding to sample
#' identifiers and columns representing summary statistics for library size
#' factors.
#'
#' The function uses the \code{\link[scater]{librarySizeFactors}} function for
#' size factor calculations and the internal function \code{.int_summaryCalc}
#' for generating individual sample summaries. This function is a crucial step
#' in understanding the distribution and characteristics of library size
#' factors in spatial transcriptomics data, providing valuable insights for
#' downstream analyses.
#'
#' @seealso \code{\link[scater]{librarySizeFactors}},
#' \code{\link{.int_summaryCalc}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics, library size factors, summary
#'
#' @rdname computeLibSizeFactors
#'
.int_summarySizeFact <- function(sfe, ids) {
  sum_list <- lapply(ids, .int_summaryCalc, sfe = sfe)
  sumOut <- Reduce(function(x, y) rbind(x, y), sum_list)
  rownames(sumOut) <- ids
  return(sumOut)
}
