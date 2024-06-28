#' Global G test for spatial autocorrelation
#'
#' This function is a wrapper of the \code{\link[spdep]{moran}} function from
#' the `spdep` package developed by Roger Bivand. The global G statistic for
#' spatial autocorrelation, complementing the local Gi LISA measures:
#' \code{\link{localG}}.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
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
#' @param alternative A character string specifying the
#'   alternative hypothesis, must be one of "greater"
#'   (default), "less" or "two.sided".
#' @param spChk Should the data vector names be checked
#'   against the spatial objects for identity integrity,
#'   TRUE, or FALSE, default NULL to use
#'   \code{get.spChkOption()}.
#' @param adjust.n Default TRUE, if FALSE the number of
#'   observations is not adjusted for no-neighbour
#'   observations, if TRUE, the number of observations is
#'   adjusted.
#' @param B1correct Default TRUE, if TRUE, the erratum
#'   referenced below: "On page 195, the coefficient of W2
#'   in B1, (just below centre of the page) should be 6,
#'   not 3." is applied; if FALSE, 3 is used (as in CrimeStat
#'   IV).
#' @param adjust.x Default TRUE, if TRUE, x values of
#'   observations with no neighbours are omitted in the
#'   denominator of G.
#' @param Arc_all_x Default FALSE, if Arc_all_x=TRUE and
#'   adjust.x=TRUE, use the full x vector in part of the
#'   denominator term for G.
#'
#' @return An SFE object with the below two columns added in the `rowData`:
#'   \item{getisG_test}{the value of the standard deviate
#'     of Getis-Ord G.}
#'   \item{getisGPval_test}{the p-value of the test.}
#'
#' @references Getis. A, Ord, J. K. 1992 The analysis of spatial association by
#' use of distance statistics, \emph{Geographical Analysis}, 24, p. 195; see
#' also Getis. A, Ord, J. K. 1993 Erratum, \emph{Geographical Analysis}, 25,
#' p. 276; Bivand RS, Wong DWS 2018 Comparing implementations of global and
#' local indicators of spatial association. TEST, 27(3), 716--748 \doi{10.1007/
#' s11749-018-0599-x}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname getisGlobalGTest
#'
#' @aliases getisGlobalGTest
#'
#' @seealso \code{\link{localG}}
#'
#' @examples
#' \dontrun{
#' # Example usage:
#'
#' # Export an SFE object from an MSFE object
#' sfe <- getSFE(msfe, "JBO019")
#'
#' # A vector of gene expression of one gene over all locations
#' gene_exp <- assay(sfe, "logcounts")[1,]
#'
#' # A listw object containing distance-based spatial weights for neighbours
#' listw <- listw
#'
#' # Calculate Moran's I for this gene
#' getis <- getisGlobalGTest(x = gene_exp, listw = listw,zero.policy = TRUE)
#' }
#'
#' @export
getisGlobalGTest <- function(m_sfe,
                             sample_id = NULL,
                             genes = TRUE,
                             zero.policy = NULL,
                             alternative = "greater",
                             spChk = NULL,
                             adjust.n = TRUE,
                             B1correct = TRUE,
                             adjust.x = TRUE,
                             Arc_all_x = FALSE,
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

  res <- parallel::mclapply(genes,
                            .int_getisCTest,
                            sfe = sfe,
                            listw = listw,
                            zero.policy = zero.policy,
                            alternative = alternative,
                            spChk = spChk,
                            adjust.n = adjust.n,
                            B1correct = B1correct,
                            adjust.x = adjust.x,
                            Arc_all_x = Arc_all_x,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's rowData
  res <- as.data.frame(rlist::list.rbind(res))
  SummarizedExperiment::rowData(sfe)$getisG_test <- unlist(res$statistic)
  SummarizedExperiment::rowData(sfe)$getisGPval_test <- unlist(res$p.value)

  ## Check and output either an msfe or an sfe object
  out <- .int_checkAndOutput(m_sfe = m_sfe, sfe = sfe, sample_id = sample_id)

  return(out)
}


#' G and Gstar local spatial statistics
#'
#' The local spatial statistic G is calculated for each
#' zone based on the spatial weights object used. The value
#' returned is a Z-value, and may be used as a diagnostic
#' tool. High positive values indicate the posibility of a
#' local cluster of high values of the variable being
#' analysed, very low relative values a similar cluster of
#' low values. For inference, a Bonferroni-type test is
#' suggested in the references, where tables of critical
#' values may be found (see also details below).
#'
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
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
#' @param spChk Should the data vector names be checked
#' against the spatial objects for identity integrity,
#' TRUE, or FALSE, default NULL to use \code{get.spChkOption()}
#' @param GeoDa Default FALSE, if TRUE, drop x values for
#' no-neighbour and self-neighbour only observations from all
#' summations
#' @param nsim Default 499, number of conditional permutation
#' simulations
#' @param alternative A character string specifying the
#' alternative hypothesis, must be one of \code{"two.sided"}
#' (default), \code{"greater"} or \code{"less"}.
#' @param iseed Default NULL, used to set the seed for
#' possible parallel RNGs
#' @param fix_i_in_Gstar_permutations Default \code{TRUE}
#' (fix x at self in permutations for local G-star), set
#' \code{FALSE} to use pre-1.2-8 behaviour
#' @param no_repeat_in_row Default \code{FALSE}, if
#' \code{TRUE}, sample conditionally in each row without
#' replacements to avoid duplicate values,
#' \url{https://github.com/r-spatial/spdep/issues/124}
#'
#' @details If the neighbours member of listw has a
#' "self.included" attribute set to TRUE, the Gstar variant,
#' including the self-weight \eqn{w_{ii} > 0}, is calculated
#' and returned. The returned vector will have a "gstari"
#' attribute set to TRUE. Self-weights must be included by
#' using the \code{include.self} function before converting
#' the neighbor list to a spatial weights list with
#' \code{nb2listw} as shown below in the example.
#'
#' The critical values of the statistic under assumptions
#' given in the references for the 95th percentile are for
#' n=1: 1.645, n=50: 3.083, n=100: 3.289, n=1000: 3.886.
#'
#' @return An SFE object with the results added in the `localResults` slot
#' which contains a DataFrame named `localGetisG`/`localGetisGPerm` that
#' contains the per location Ci statistics values for each gene:
#' \item{ENSG***.Gi}{A vector of G or Gstar statistic (standard deviate values)}
#' \item{ENSG***.GiFDR}{P-value of the calculated statistic}
#' \item{ENSG***.GiClust}{High positive values indicate the posibility of a
#' local cluster of high values of the variable being
#' analysed, very low relative values a similar cluster of
#' low values.}
#'
#' @note Conditional permutations added for comparative
#' purposes; permutations are over the whole data vector
#' omitting the observation itself, and from 1.2-8 fixing the
#' observation itself as its own neighbour for local G-star.
#'
#' @references Ord, J. K. and Getis, A. 1995 Local spatial
#' autocorrelation statistics: distributional issues and an
#' application. \emph{Geographical Analysis}, 27, 286--306;
#' Getis, A. and Ord, J. K. 1996 Local spatial statistics: an
#' overview. In P. Longley and M. Batty (eds) \emph{Spatial
#' analysis: modelling in a GIS environment} (Cambridge:
#' Geoinformation International), 261--277; Bivand RS, Wong
#' DWS 2018 Comparing implementations of global and local
#' indicators of spatial association.
#' TEST, 27(3), 716--748 \doi{10.1007/s11749-018-0599-x}
#'
#' @seealso \code{\link{localG}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname getisLocalG
#'
#' @aliases getisLocalG
#'
#' @examples
#' \dontrun{
#' # Example usage:
#'
#' # Export an SFE object from an MSFE object
#' sfe <- getSFE(msfe, "JBO019")
#'
#' # A vector of gene expression of one gene over all locations
#' gene_exp <- assay(sfe, "logcounts")[1,]
#'
#' # A listw object containing distance-based spatial weights for neighbours
#' listw <- listw
#'
#' # Calculate Getis and Ord's G for this gene
#' getis <- getisLocalG(x = gene_exp, listw = listw,zero.policy = TRUE)
#' }
#'
#' @export
getisLocalG <- function(m_sfe,
                        sample_id = NULL,
                        genes = TRUE,
                        zero.policy = NULL,
                        spChk = NULL,
                        GeoDa = FALSE,
                        alternative = "two.sided",
                        return_internals = TRUE,
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

  res <- parallel::mclapply(genes,
                            .int_getisLocal,
                            sfe = sfe,
                            listw = listw,
                            zero.policy = zero.policy,
                            spChk = spChk,
                            GeoDa = GeoDa,
                            alternative = alternative,
                            return_internals = return_internals,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's localResults
  res <- S4Vectors::DataFrame(rlist::list.cbind(res))
  localResults(sfe, name = "localGetisG") <- res

  ## Check and output either an msfe or an sfe object
  out <- .int_checkAndOutput(m_sfe = m_sfe, sfe = sfe, sample_id = sample_id)

  return(out)
}



#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname getisLocalG
#'
#' @aliases getisLocalGPerm
#'
#' @seealso \code{\link{localG}}
#'
#' @examples
#' \dontrun{
#' # Calculate Getis and Ord's G for this gene
#' getis <- getisLocalGPerm(x = gene_exp, listw = listw, nsim = 499,
#' zero.policy = TRUE)
#' }
#'
#' @export
getisLocalGPerm <- function(m_sfe,
                            sample_id = NULL,
                            genes = TRUE,
                            nsim = 999,
                            zero.policy = NULL,
                            spChk = NULL,
                            alternative = "two.sided",
                            iseed = NULL,
                            fix_i_in_Gstar_permutations = TRUE,
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
    stop("Invalid `genes` argument input")
  }

  ## Get neighbour graph
  listw <- colGraph(sfe)

  ## Check zero.policy
  if (is.null(zero.policy)) {
    zero.policy = attr(listw, "zero.policy")
  }

  res <- parallel::mclapply(genes,
                            .int_getisLocalPerm,
                            sfe = sfe,
                            listw = listw,
                            nsim = nsim,
                            zero.policy = zero.policy,
                            spChk = spChk,
                            alternative = alternative,
                            iseed = iseed,
                            fix_i_in_Gstar_permutations =
                              fix_i_in_Gstar_permutations,
                            no_repeat_in_row = no_repeat_in_row,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's localResults
  res <- S4Vectors::DataFrame(rlist::list.cbind(res))
  localResults(sfe, name = "localGetisGPerm") <- res

  ## Check and output either an msfe or an sfe object
  out <- .int_checkAndOutput(m_sfe = m_sfe, sfe = sfe, sample_id = sample_id)

  return(out)
}


# ---------------------------------------------------------------------------- #
#  ###### INTERNAL FUNCTIONS ASSOCIATED WITH SA CALCULATIONS (C, G, I) ######
# ---------------------------------------------------------------------------- #
#' Internal Function: Calculate Getis & Ord's Global G with statistical
#' significance
#'
#' This internal function calculates the G statistic and the statistical
#' significance for a vector of gene expression. For the parameter arguments
#' check the \code{link[STExplorer]{getisGlobalGTest}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getisGTest
#'
#' @aliases .int_getisGTest
#'
.int_getisGTest <- function(gene,
                            sfe,
                            listw,
                            zero.policy,
                            alternative,
                            spChk,
                            adjust.n,
                            B1correct,
                            adjust.x,
                            Arc_all_x) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the globalG.test function from spdep
  spdep::globalG.test(x = x,
                      listw = listw,
                      zero.policy = zero.policy,
                      alternative = alternative,
                      spChk = spChk,
                      adjust.n = adjust.n,
                      B1correct = B1correct,
                      adjust.x = adjust.x,
                      Arc_all_x = Arc_all_x)
}


#' Internal Function: Calculate Getis & Ord's G
#'
#' This internal function calculates the G statistic for a vector of gene
#' expression. For the parameter arguments check the
#' \code{link[STExplorer]{getisGlobalG}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getisLocal
#'
#' @aliases .int_getisLocal
#'
.int_getisLocal <- function(gene,
                            sfe,
                            listw,
                            zero.policy,
                            spChk,
                            GeoDa,
                            alternative,
                            return_internals) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localG function from spdep
  res <- spdep::localG(x = x,
                       listw = listw,
                       zero.policy = zero.policy,
                       spChk = spChk,
                       GeoDa = GeoDa,
                       alternative = alternative,
                       return_internals = return_internals)

  ## Get attributes into a dataframe
  fdrColName <- grep("Pr.*",
                     colnames(attr(res, "internals")),
                     value = TRUE)
  res <- data.frame(Gi = as.vector(res),
                    GiFDR = attr(res, "internals")[,fdrColName],
                    GiClust = attr(res, "cluster"))

  return(res)
}


#' Internal Function: Calculate Geary's C
#'
#' This internal function calculates the C statistic for a vector of gene
#' expression. For the parameter arguments check the
#' \code{link[STExplorer]{gearyGlobalC}} function.

#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getisLocal
#'
#' @aliases .int_getisLocal
#'
.int_getisLocalPerm <- function(gene,
                                sfe,
                                listw,
                                nsim,
                                zero.policy,
                                spChk,
                                alternative,
                                iseed,
                                fix_i_in_Gstar_permutations,
                                no_repeat_in_row) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localG_perm function from spdep
  res <- spdep::localG_perm(x = x,
                            listw = listw,
                            nsim = nsim,
                            zero.policy = zero.policy,
                            spChk = spChk,
                            alternative = alternative,
                            iseed = iseed,
                            fix_i_in_Gstar_permutations = fix_i_in_Gstar_permutations,
                            no_repeat_in_row = no_repeat_in_row)

  ## Get attributes into a dataframe
  fdrColName <- grep("Pr.z.*Sim",
                     colnames(attr(res, "internals")),
                     value = TRUE)
  res <- data.frame(Gi = as.vector(res),
                    GiFDR = attr(res, "internals")[,fdrColName],
                    GiClust = attr(res, "cluster"))

  return(res)
}
