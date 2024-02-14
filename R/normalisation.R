#' Compute Library Size Factors in a SpatialFeatureExperiment
#'
#' This function computes library size factors in a
#' \code{SpatialFeatureExperiment} object. Size factors represent a scaling
#' factor for each sample, adjusting for differences in sequencing depth.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character string specifying unique sample identifiers,
#' one for each directory specified via samples. This parameter is ignored if
#' the names of samples are not null (!is.null(names(samples))).
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
computeLibSizeFactors <- function(m_sfe,
                                  sample_id,
                                  type = c("inter", "intra"),
                                  ...) {
  UseMethod("computeLibSizeFactors")
}


#' @export
computeLibSizeFactors.SpatialFeatureExperiment <-
  function(m_sfe,
           sample_id,
           type = c("inter", "intra"),
           ...) {

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
    m_sfe <- scater::computeLibraryFactors(m_sfe, ...)

    ## Return the value
    return(m_sfe)
  } else if (in_sample) {
    ## Calculate library size factors
    ## Take into account one sample at a time

    ## Get sample IDs
    ids <- .int_getSmplIDs(sfe = m_sfe, sample_id = TRUE)

    ## Create a list of filtered SpatialFeatureExperiment objects
    sfe_list <- lapply(ids, .int_sizeFactCalc, sfe = m_sfe, ...)

    ## Merge the list of SpatialFeatureExperiment objects into one
    sfeOut <- Reduce(function(x, y) cbind(x, y), sfe_list)

    ## Clean up duplicate entries in the metadata slot
    sfeOut <- .int_cleanMetaData(sfeOut = sfeOut)

    ## Have a look at the size factors
    print(.int_summarySizeFact(sfe = sfeOut, ids = ids))

    ## Return the value
    return(sfeOut)
  }
}


#' @export
computeLibSizeFactors.MetaSpatialFeatureExperiment <-
  function(m_sfe,
           sample_id = TRUE,
           type = NULL,
           ...) {

  ## Check valid type argument
  if (is.null(type)) {
    in_sample <- TRUE
  } else {
    type <- match.arg(type)
    if (type == "intra") {
      in_sample <- TRUE
    } else {
      in_sample <- FALSE
    }
  }

  if (sample_id) {
    ## Calculate library size factors
    ## Take into account all samples together
    m_sfe@sfe_data <- lapply(m_sfe@sfe_data,
                             .int_msfeCompSizeFact,
                             ... = ...)

    return(m_sfe)
  } else if (is.character(sample_id)) {
    ## Calculate library size factors
    ## Take into account only specified samples
    m_sfe@sfe_data <- lapply(m_sfe@sfe_data[[sample_id]],
                             .int_msfeCompSizeFact,
                             ... = ...)

    return(m_sfe)
  }
}


#' Normalise Counts in a Spatial Transcriptomics Experiment
#'
#' This function normalises counts in a spatial transcriptomics experiment
#' using specified methods.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character string specifying unique sample identifiers,
#' one for each directory specified via samples. This parameter is ignored if
#' the names of samples are not null (!is.null(names(samples))).
#' @param method The normalisation method to be applied. Currently, only "log"
#' transformation is supported.
#' @param ... Additional arguments to be passed to the normalisation function.
#'
#' @return A normalised SpatialFeatureExperiment object.
#'
#' @details The function normalises counts in a spatial transcriptomics
#' experiment using specified methods. Currently, only the "log" transformation
#' method is supported. Additional arguments can be passed to the underlying
#' normalisation function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords normalisation SpatialFeatureExperiment log transformation
#'
#' @rdname normaliseCounts
#'
#' @examples
#' \dontrun{
#' # Normalise counts using log transformation
#' # msfe <- # MetaSpatialFeatureExperiment object
#' msfe <- normaliseCounts(msfe, method = "log", base = 2)
#' }
normaliseCounts <- function(m_sfe,
                            sample_id,
                            method = c("log"),
                            ...) {
  UseMethod("normaliseCounts")
}


#' @export
normaliseCounts.SpatialFeatureExperiment <- function(m_sfe,
                                                     method = c("log"),
                                                     ...) {
  ## Check valid method argument
  if (missing(method)) {
    method <- "log"
  } else {
    method <- match.arg(method)
  }

  if (method == "log") {
    m_sfe <- lapply(m_sfe, scater::logNormCounts, ...)
  }

  return(m_sfe)
}


#' @export
normaliseCounts.MetaSpatialFeatureExperiment <- function(m_sfe,
                                                         method = c("log"),
                                                         ...) {
  ## Check valid method argument
  if (missing(method)) {
    method <- "log"
  } else {
    method <- match.arg(method)
  }

  if (method == "log") {
    m_sfe@sfe_data <- lapply(msfe@sfe_data, scater::logNormCounts, ...)
  }

  return(m_sfe)
}


# ---------------------------------------------------------------------------- #
#  ########## INTERNAL FUNCTIONS ASSOCIATED WITH LIB SIZE FACTORS ###########
# ---------------------------------------------------------------------------- #
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


# ---------------------------------------------------------------------------- #
#  ########## INTERNAL FUNCTIONS ASSOCIATED WITH LIB SIZE FACTORS ###########
# ---------------------------------------------------------------------------- #
#' INTERNAL: Compute Library Size Factors in a SpatialFeatureExperiment
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param ... other arguments to be passed to `scater::computeLibraryFactors`
#'
#' @return A logical vector indicating the subset of features based on the
#' criteria.
#'
#' @rdname dot-int_msfeCompSizeFact
#'
.int_msfeCompSizeFact <- function(sfe, ...) {
  scater::computeLibraryFactors(sfe, ...)
}
