#' Calculate Geary's C Global Statistic
#'
#' This function is a wrapper of the \code{\link[spdep]{geary}} function from
#' the `spdep` package developed by Roger Bivand. It computes Geary's C
#' statistic for spatial autocorrelation using the spdep package. In its core
#' the function asks for two arguments: 1) `x` - A numeric vector of
#' observations. Must be the same length as the neighbours list in listw, and
#' 2) `listw` A listw object created, for example, by
#' \code{link[spdep]{nb2listw}} or \code{link[spdep]{nb2listwdist}}.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
#' @param n Number of zones. If not provided, it is calculated as
#' length(listw$neighbours).
#' @param n1 \code{n - 1}. If not provided, it is calculated as
#' \code{length(listw$neighbours) - 1}.
#' @param S0 Global sum of weights. If not provided, it is calculated as
#' \code{spdep::Szero(listw)}.
#' @param zero.policy Default is NULL. If not changed then internally, it is
#' set to \code{attr(listw, "zero.policy")} as set when \code{listw} was
#' created. If TRUE, assign zero to the lagged value of zones without
#' neighbours. If FALSE, assign NA. If attribute not set, use the global option
#' value.
#' @param mc.cores Argument from \code{link[parallel]{mclapply}}. The number of
#' cores to use, i.e., at most how many child processes will be run
#' simultaneously. The option is initialized from environment variable MC_CORES
#' if set. Must be at least one, and parallelisation requires at least two
#' cores.
#'
#' @return An SFE object with the below two columns added in the `rowData`:
#' \itemize{
#'  \item{gearyC}{the value of the observed Geary's C.}
#'  \item{gearyK}{sample kurtosis of gene.}
#' }
#'
#' @details If n, n1, and S0 are not provided, they are calculated based on the
#' \code{listw} object.
#'
#' @seealso \code{\link[spdep]{geary}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords geary spatial autocorrelation listw
#'
#' @rdname gearyGlobalC
#'
#' @aliases gearyGlobalC
#'
#' @examples
#' \dontrun{
#' # Export an SFE object from an MSFE object
#' sfe <- getSFE(msfe, "JBO019")
#' # A vector of gene expression of one gene over all locations
#' gene_exp <- assay(sfe, "logcounts")[1,]
#' # A listw object containing distance-based spatial weights for neighbors
#' listw <- listw
#' # Calculate Geary's C for this gene
#' geary <- gearyGlobalC(x = gene_exp, listw = listw, zero.policy = TRUE)
#' }
#'
#' @export
gearyGlobalC <- function(m_sfe,
                         sample_id,
                         genes = TRUE,
                         n = NULL,
                         n1 = NULL,
                         S0 = NULL,
                         zero.policy = NULL,
                         mc.cores = getOption("mc.cores", 2L)) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check genes to use
  if (is.character(genes)) {
    # genes is already a character vector, no need to modify it
  } else if (isTRUE(genes)) {
    genes <- rownames(rowData(sfe))
    names(genes) <- genes
  } else {
    stop("Invalid `genes` argument input")
  }

  ## Get neighbour graph
  listw <- colGraph(sfe)

  ## Check zero.policy
  if (is.null(zero.policy)) {
    zero.policy = attr(listw, "zero.policy")
  }

  ## Calculate n, n1, and S0 if not provided
  if (is.null(n)) n <- length(listw$neighbours)
  if (is.null(n1)) n1 <- n - 1
  if (is.null(S0)) S0 <- spdep::Szero(listw)

  ## Call the geary function from spdep
  out <- parallel::mclapply(genes,
                  .int_geary,
                  sfe = sfe,
                  listw = listw,
                  n = n,
                  n1 = n1,
                  S0 = S0,
                  zero.policy = zero.policy,
                  mc.cores = mc.cores)

  ## Import output into the SFE object's rowData
  out <- as.data.frame(rlist::list.rbind(out))
  SummarizedExperiment::rowData(sfe)$gearyC <- out$C
  SummarizedExperiment::rowData(sfe)$gearyK <- out$K

  return(sfe)
}

#' Calculate Geary's C Global Statistic with permutation testing
#'
#' This function is a wrapper of the \code{\link[spdep]{geary.mc}} function
#' from the `spdep` package developed by Roger Bivand. It computes Geary's C
#' statistic with permutation testing for spatial autocorrelation using the
#' spdep package. A permutation test for Geary's C statistic calculated by
#' using `nsim` random permutations of `x` for the given spatial weighting
#' scheme, to establish the rank of the observed statistic in relation to the
#' `nsim` simulated values.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
#' @param nsim Number of permutations.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "greater" (default), or "less"; this reversal corresponds to
#' that on geary.test described in the section on the output statistic value,
#' based on Cliff and Ord 1973, p. 21.
#' @param spChk should the data vector names be checked against the spatial
#' objects for identity integrity, TRUE, or FALSE, default NULL to use
#' \code{get.spChkOption()}.
#' @param zero.policy Default is NULL. If not changed then internally, it is
#' set to \code{attr(listw, "zero.policy")} as set when \code{listw} was
#' created. If TRUE, assign zero to the lagged value of zones without
#' neighbours. If FALSE, assign NA. If attribute not set, use the global option
#' value.
#' @param adjust.n Default TRUE, if FALSE the number of observations is not
#' adjusted for no-neighbour observations, if TRUE, the number of observations
#' is adjusted.
#' @param return_boot Return an object of class boot from the equivalent
#' permutation bootstrap rather than an object of class htest.
#' @param mc.cores Argument from \code{link[parallel]{mclapply}}. The number of
#' cores to use, i.e., at most how many child processes will be run
#' simultaneously. The option is initialized from environment variable MC_CORES
#' if set. Must be at least one, and parallelisation requires at least two
#' cores.
#'
#'
#' @return An SFE object with the below two columns added in the `rowData`:
#' \itemize{
#'  \item{gearyC_perm}{the value of the observed Geary's C.}
#'  \item{gearyPval_perm}{the pseudo p-value of the test.}
#' }
#'
#' @details If n, n1, and S0 are not provided, they are calculated based on the
#' \code{listw} object.
#'
#' @seealso \code{\link[spdep]{geary.test}} \code{\link[spdep]{geary}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords geary spatial autocorrelation listw
#'
#' @rdname gearyGlobalCPerm
#'
#' @aliases gearyGlobalCPerm
#'
#' @examples
#' \dontrun{
#' # Export an SFE object from an MSFE object
#' sfe <- getSFE(msfe, "JBO019")
#' # A vector of gene expression of one gene over all locations
#' gene_exp <- assay(sfe, "logcounts")[1,]
#' # A listw object containing distance-based spatial weights for neighbours
#' listw <- listw
#' # Calculate Geary's C and statistical significance for this gene
#' geary <- gearyGlobalCPerm(x = gene_exp, listw = listw, nsim = 999,
#' zero.policy = TRUE)
#' }
#'
#' @export
gearyGlobalCPerm <- function(m_sfe,
                             sample_id,
                             genes = TRUE,
                             nsim,
                             zero.policy = NULL,
                             alternative = "greater",
                             spChk = NULL,
                             adjust.n = TRUE,
                             return_boot = FALSE,
                             mc.cores = getOption("mc.cores", 2L)) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check genes to use
  if (is.character(genes)) {
    # genes is already a character vector, no need to modify it
  } else if (isTRUE(genes)) {
    genes <- rownames(rowData(sfe))
    names(genes) <- genes
  } else {
    stop("Invalid `genes` argument input")
  }

  ## Get neighbour graph
  listw <- colGraph(sfe)

  ## Assign nsim if not provided
  if (missing(nsim)) {
    nsim <- 999
  }

  ## Check zero.policy
  if (is.null(zero.policy)) {
    zero.policy = attr(listw, "zero.policy")
  }

  ## Call the geary function from spdep
  out <- parallel::mclapply(genes,
                            .int_gearyCPerm,
                            sfe = sfe,
                            listw = listw,
                            nsim = nsim,
                            alternative = alternative,
                            spChk = spChk,
                            adjust.n = adjust.n,
                            return_boot = return_boot,
                            zero.policy = zero.policy,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's rowData
  out <- as.data.frame(rlist::list.rbind(out))
  SummarizedExperiment::rowData(sfe)$gearyC_perm <- out$statistic
  SummarizedExperiment::rowData(sfe)$gearyPval_perm <- out$p.value

  return(out)
}

#' Geary's C Test for Spatial Autocorrelation
#'
#' This function is a wrapper of the \code{\link[spdep]{geary.test}} function
#' from the `spdep` package developed by Roger Bivand. It performs Geary's test
#' for spatial autocorrelation using a spatial weights matrix in weights list
#' form. The assumptions underlying the test are sensitive to the form of the
#' graph of neighbour relationships and other factors, and results may be
#' checked against those of \code{\link[spdep]{geary.mc}} permutations.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
#' @param randomisation Variance of I calculated under the assumption of
#' randomisation. If FALSE, normality.
#' @param zero.policy Default is NULL. If not changed then internally, it is
#' set to \code{attr(listw, "zero.policy")} as set when \code{listw} was
#' created. If TRUE, assign zero to the lagged value of zones without
#' neighbours. If FALSE, assign NA. If attribute not set, use the global option
#' value.
#' @param alternative A character string specifying the alternative hypothesis.
#' Must be one of "greater" (default), "less" or "two.sided".
#' @param spChk Should the data vector names be checked against the spatial
#' objects for identity integrity? TRUE or FALSE, default NULL to use
#' \code{get.spChkOption()}.
#' @param adjust.n Default TRUE. If FALSE, the number of observations is not
#' adjusted for no-neighbour observations. If TRUE, the number of observations
#' is adjusted.
#' @param mc.cores Argument from \code{link[parallel]{mclapply}}. The number of
#' cores to use, i.e., at most how many child processes will be run
#' simultaneously. The option is initialized from environment variable MC_CORES
#' if set. Must be at least one, and parallelisation requires at least two
#' cores.
#'
#' @return An SFE object with the below two columns added in the `rowData`:
#'   \item{gearyC_test}{the value of the standard deviate of Geary's C.}
#'   \item{gearyCPval_test}{the p-value of the test.}
#'
#' @details The derivation of the test (Cliff and Ord, 1981, p. 18) assumes
#' that the weights matrix is symmetric. For inherently non-symmetric matrices,
#' such as k-nearest neighbour matrices, \code{listw2U()} can be used to make
#' the matrix symmetric. In non-symmetric weights matrix cases, the variance of
#' the test statistic may be negative. Geary's C is affected by non-symmetric
#' weights under normality much more than Moran's I. From 0.4-35, the sign of
#' the standard deviate of C is changed to match Cliff and Ord (1973, p. 21).
#'
#' @references
#' Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21,
#' Cliff, A. D., Ord, J. K. 1973 Spatial Autocorrelation, Pion, pp. 15-16, 21;
#' Bivand RS, Wong DWS 2018 Comparing implementations of global and local
#' indicators of spatial association. TEST, 27(3), 716--748 \doi{10.1007/s11749
#' -018-0599-x}
#'
#' @seealso \code{\link{geary}}, \code{\link{geary.mc}}, \code{\link{listw2U}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords geary spatial autocorrelation listw
#'
#' @rdname gearyGlobalCTest
#'
#' @aliases gearyGlobalCTest
#'
#' @examples
#' \dontrun{
#' # Export an SFE object from an MSFE object
#' sfe <- getSFE(msfe, "JBO019")
#' # A vector of gene expression of one gene over all locations
#' gene_exp <- assay(sfe, "logcounts")[1,]
#' # A listw object containing distance-based spatial weights for neighbours
#' listw <- listw
#' # Calculate Geary's C and statistical significance for this gene
#' geary <- gearyGlobalCTest(x = gene_exp, listw = listw, zero.policy = TRUE)
#' }
#'
#' @export
gearyGlobalCTest <- function(m_sfe,
                             sample_id,
                             genes = TRUE,
                             randomisation=TRUE,
                             zero.policy = NULL,
                             alternative = "greater",
                             spChk = NULL,
                             adjust.n = TRUE,
                             mc.cores = getOption("mc.cores", 2L)) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check genes to use
  if (is.character(genes)) {
    # genes is already a character vector, no need to modify it
  } else if (isTRUE(genes)) {
    genes <- rownames(rowData(sfe))
    names(genes) <- genes
  } else {
    stop("Invalid `genes` argument input")
  }

  ## Get neighbour graph
  listw <- colGraph(sfe)

  ## Check zero.policy
  if (is.null(zero.policy)) {
    zero.policy = attr(listw, "zero.policy")
  }

  ## Call the geary.test function from spdep
  out <- parallel::mclapply(genes,
                            .int_gearyTest,
                            sfe = sfe,
                            listw = listw,
                            randomisation = randomisation,
                            alternative = alternative,
                            spChk = spChk,
                            adjust.n = adjust.n,
                            zero.policy = zero.policy,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's rowData
  out <- as.data.frame(rlist::list.rbind(out))
  SummarizedExperiment::rowData(sfe)$gearyC_test <- out$statistic
  SummarizedExperiment::rowData(sfe)$gearyCPval_test <- out$p.value

  return(sfe)
}


#' Local Geary's C Statistic for Spatial Autocorrelation
#'
#' The Local Geary is a local adaptation of Geary's C
#' statistic of spatial autocorrelation. It uses squared
#' differences to measure dissimilarity, unlike the Local
#' Moran. Low values indicate positive spatial autocorrelation,
#' and large values refer to negative autocorrelation.
#'
#' Inference for the Local Geary is based on permutation
#' approach, comparing the observed value to the reference
#' distribution under spatial randomness. \code{localC_perm()}
#' returns a pseudo p-value. This is not an analytical p-value
#' and should be used with care, based on permutations.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
#' @param nsim The number of simulations to be used for permutation test.
#' @param alternative A character defining the alternative hypothesis. Must be
#' one of "two.sided", "less" or "greater".
#' @param ... other arguments passed to methods. Have a look at
#' \code{\link[spdep]{localC}}
#' @param iseed default NULL, used to set the seed for possible parallel RNGs.
#' @param no_repeat_in_row default FALSE, if TRUE, sample conditionally in each
#' row without replacements to avoid duplicate values,
#' https://github.com/r-spatial/spdep/issues/124
#' @param zero.policy Default is NULL. If not changed then internally, it is
#' set to \code{attr(listw, "zero.policy")} as set when \code{listw} was
#' created. If TRUE, assign zero to the lagged value of zones without
#' neighbours. If FALSE, assign NA. If attribute not set, use the global option
#' value.
#' @param mc.cores Argument from \code{link[parallel]{mclapply}}. The number of
#' cores to use, i.e., at most how many child processes will be run
#' simultaneously. The option is initialized from environment variable MC_CORES
#' if set. Must be at least one, and parallelisation requires at least two
#' cores.
#'
#' @details For more information about the usage cases and the
#' rest of the arguments visit the help page of \code{\link[spdep]{localC}}
#' Local Geary can be extended to a multivariate
#' context. When \code{x} is numeric, univariate Local Geary
#' is calculated. For multivariate, provide a list or matrix.
#'
#' While not required in univariate context, standardized
#' Local Geary is calculated. Multivariate Local Geary is
#' \emph{always} standardized.
#'
#' The univariate Local Geary is calculated as \eqn{c_i = \sum_j
#' w_{ij}(x_i - x_j)^2} and the multivariate Local Geary is
#' calculated as \eqn{c_{k,i} = \sum_{v=1}^{k} c_{v,i}} as
#' described in Anselin (2019).
#'
#' @return An SFE object with the results added in the `localResults` slot
#' which contains a DataFrame named `localGearyC`/`localGearyCPerm` that
#' contains the per location Ci statistics values for each gene:
#' \item{Ci or ENSG***.Ci}{Numeric vector containing Local Geary statistic
#' with attribute
#' \item{ENSG***.CiFDR}{(Only for gearyLocalCPerm) \code{rank()} and
#' \code{punif()} of observed statistic rank for [0, 1] p-values using
#' \code{alternative=}}
#' \item{ENSG***.CiCluster}{(Only for gearyLocalCPerm) Low/High cluster for
#' each location. It can take values: "Low-Low", "High-High", "Other Positive",
#' "Negative"}
#'
#' @references Anselin, L. (2019). Local indicators of spatial
#' association—LISA. Spatial Statistics, 29, 10-27.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords geary spatial autocorrelation listw
#'
#' @rdname gearyLocalC
#'
#' @aliases gearyLocalC
#'
#' @examples
#' \dontrun{
#' # Export an SFE object from an MSFE object
#' sfe <- getSFE(msfe, "JBO019")
#' # A vector of gene expression of one gene over all locations
#' gene_exp <- assay(sfe, "logcounts")[1,]
#' # A listw object containing distance-based spatial weights for neighbours
#' listw <- listw
#' # Calculate Geary's C for this gene
#' geary <- gearyLocalC(x = gene_exp, listw = listw, zero.policy = TRUE)
#' }
#'
#' @export
gearyLocalC <- function(m_sfe,
                        sample_id,
                        genes = TRUE,
                        ...,
                        zero.policy = NULL,
                        mc.cores = getOption("mc.cores", 2L)) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check genes to use
  if (is.character(genes)) {
    # genes is already a character vector, no need to modify it
  } else if (isTRUE(genes)) {
    genes <- rownames(rowData(sfe))
    names(genes) <- genes
  } else {
    stop("Invalid `genes` argument input")
  }

  ## Get neighbour graph
  listw <- colGraph(sfe)

  ## Check zero.policy
  if (is.null(zero.policy)) {
    zero.policy = attr(listw, "zero.policy")
  }

  out <- parallel::mclapply(genes,
                            .int_gearyLocalC,
                            sfe = sfe,
                            listw = listw,
                            ...,
                            zero.policy = zero.policy,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's localResults
  out <- S4Vectors::DataFrame(rlist::list.cbind(out))
  localResults(sfe, name = "localGearyC") <- out

  return(sfe)
}


#' @aliases gearyLocalCPerm
#' @rdname gearyLocalC
#' @return description
#' @examples
#' \dontrun{
#' # Calculate  Geary's C and statistical significance for this gene
#' geary <- gearyLocalCPerm(x = gene_exp, listw = listw, nsim = 499,
#' zero.policy = TRUE)
#' }
#' @export
gearyLocalCPerm <- function(m_sfe,
                            sample_id,
                            genes = TRUE,
                            nsim = 999,
                            alternative = "two.sided",
                            ...,
                            zero.policy = NULL,
                            iseed = NULL,
                            no_repeat_in_row = FALSE,
                            mc.cores = getOption("mc.cores", 2L)) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check genes to use
  if (is.character(genes)) {
    # genes is already a character vector, no need to modify it
  } else if (isTRUE(genes)) {
    genes <- rownames(rowData(sfe))
    names(genes) <- genes
  } else {
    stop("Invalid `genes` argument input.")
  }

  ## Get neighbour graph
  listw <- colGraph(sfe)

  ## Check zero.policy
  if (is.null(zero.policy)) {
    zero.policy = attr(listw, "zero.policy")
  }

  out <- parallel::mclapply(genes,
                            .int_gearyLocalCPerm,
                            sfe = sfe,
                            listw = listw,
                            nsim = nsim,
                            alternative = alternative,
                            ...,
                            zero.policy = zero.policy,
                            iseed = iseed,
                            no_repeat_in_row = no_repeat_in_row,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's localResults
  out <- S4Vectors::DataFrame(rlist::list.cbind(out))
  localResults(sfe, name = "localGearyCPerm") <- out

  return(sfe)
}


# ---------------------------------------------------------------------------- #
#  ###### INTERNAL FUNCTIONS ASSOCIATED WITH SA CALCULATIONS (C, G, I) ######
# ---------------------------------------------------------------------------- #
#' Internal Function: Check Spatial Autocorrelation Input
#'
#' This internal function checks the validity of spatial autocorrelation input.
#'
#' @param x A numeric vector, the same length as the neighbours list in listw.
#' @param listw A \code{listw} object created, for example, by \code{nb2listw}.
#'
#' @details
#' This function is intended for internal use to validate the input parameters
#' for spatial autocorrelation calculations. It checks if 'x' is numeric, if
#' its length matches the length of neighbors in 'listw', and if 'listw' is of
#' class 'listw'.
#'
#' @return
#' This function does not return any value. It stops and throws an error if the
#' input is invalid.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_checkSAInput
#'
#' @aliases .int_checkSAInput
.int_checkSAInput <- function(x, listw) {
  if (!is.numeric(x) ||
      length(x) != length(listw$neighbours) ||
      !inherits(listw, "listw")) {
    stop("Invalid input. Please provide a numeric vector 'x', listw",
         " of class 'listw', and ensure the length of 'x' matches the",
         " length of listw$neighbours.")
  }
}


#' Internal Function: Calculate Geary's C
#'
#' This internal function calculates the C statistic for a vector of gene
#' expression. For the parameter arguments check the
#' \code{link[STExplorer]{gearyGlobalC}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @importFrom SummarizedExperiment assay
#'
#' @rdname dot-int_geary
#'
#' @aliases .int_geary
.int_geary <- function(gene,
                       sfe,
                       listw,
                       n,
                       n1,
                       S0,
                       zero.policy) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the geary function from spdep
  spdep::geary(x = x,
               listw = listw,
               n = n,
               n1 = n1,
               S0 = S0,
               zero.policy = zero.policy)
}


#' Internal Function: Calculate Geary's C with Permutation testing
#'
#' This internal function calculates the C statistic for a vector of gene
#' expression with Monte-Carlo permutation testing of significance. For the
#' parameter arguments check the \code{link[STExplorer]{gearyGlobalCPerm}}
#' function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_gearyCPerm
#'
#' @aliases .int_gearyCPerm
#'
.int_gearyCPerm <- function(gene,
                         sfe,
                         listw,
                         nsim,
                         alternative,
                         spChk,
                         adjust.n,
                         return_boot,
                         zero.policy) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the geary function from spdep
  spdep::geary.mc(x = x,
                  listw = listw,
                  nsim = nsim,
                  alternative = alternative,
                  spChk = spChk,
                  adjust.n = adjust.n,
                  return_boot = return_boot,
                  zero.policy = zero.policy)
}


#' Internal Function: Calculate Geary's C with statistical testing
#'
#' This internal function calculates the C statistic for a vector of gene
#' expression alongside a statistical significance testing. For the parameter
#' arguments check the \code{link[STExplorer]{gearyGlobalCTest}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_gearyCTest
#'
#' @aliases .int_gearyCTest
#'
.int_gearyCTest <- function(gene,
                            sfe,
                            listw,
                            randomisation,
                            alternative,
                            spChk,
                            adjust.n,
                            zero.policy) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the geary function from spdep
  spdep::geary.test(x = x,
                    listw = listw,
                    randomisation = randomisation,
                    alternative = alternative,
                    spChk = spChk,
                    adjust.n = adjust.n,
                    zero.policy = zero.policy)
}


#' Internal Function: Calculate local Geary's C
#'
#' This internal function calculates the local Ci statistic for a vector of
#' gene expression. For the parameter arguments check the
#' \code{link[STExplorer]{gearyLocalC}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_gearyLocalC
#'
#' @aliases .int_gearyLocalC
#'
.int_gearyLocalC <- function(gene,
                             sfe,
                             listw,
                             zero.policy,
                             ...) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

 res <-  spdep::localC(x = x,
                listw = listw,
                ...,
                zero.policy = zero.policy)

  out <- data.frame(Ci = as.vector(res))

  return(out)
}


#' Internal Function: Calculate local Geary's C with statistical testing
#'
#' This internal function calculates the local Ci statistic for a vector of
#' gene expression alongside a statistical significance testing. For the
#' parameter arguments check the \code{link[STExplorer]{gearyLocalCPerm}}
#' function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_gearyLocalCPerm
#'
#' @aliases .int_gearyLocalCPerm
#'
.int_gearyLocalCPerm <- function(gene,
                                 sfe,
                                 listw,
                                 nsim,
                                 alternative,
                                 ...,
                                 zero.policy,
                                 iseed,
                                 no_repeat_in_row) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  res <- spdep::localC_perm(x = x,
                            listw = listw,
                            nsim = nsim,
                            alternative = alternative,
                            ...,
                            zero.policy = zero.policy,
                            iseed = iseed,
                            no_repeat_in_row = no_repeat_in_row)

  ## Get attributes into a dataframe
  fdrColName <- grep("Pr.z.*Sim",
                     colnames(attr(res, "pseudo-p")),
                     value = TRUE)
  out <- data.frame(Ci = as.vector(res),
                    CiFDR = attr(res, "pseudo-p")[,fdrColName],
                    CiClust = attr(res, "cluster"))

  return(out)
}
