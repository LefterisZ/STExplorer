#' Global G test for spatial autocorrelation
#'
#' This function is a wrapper of the \code{\link[spdep]{moran}} function from
#' the `spdep` package developed by Roger Bivand. The global G statistic for
#' spatial autocorrelation, complementing the local Gi LISA measures:
#' \code{\link{localG}}.
#'
#' @param x Numeric vector the same length as the neighbors
#'   list in listw
#' @param listw \code{listw} object created, for example,
#'   by \code{nb2listw}; if a sequence of distance bands is
#'   to be used, it is recommended that the weights style be
#'   binary (one of \code{c("B", "C", "U")}).
#' @param zero.policy Default \code{attr(listw,
#'   "zero.policy")} as set when \code{listw} was created. If
#'   attribute not set, use global option value. If TRUE,
#'   assign zero to the lagged value of zones without
#'   neighbors; if FALSE, assign NA.
#' @param alternative A character string specifying the
#'   alternative hypothesis, must be one of "greater"
#'   (default), "less" or "two.sided".
#' @param spChk Should the data vector names be checked
#'   against the spatial objects for identity integrity,
#'   TRUE, or FALSE, default NULL to use
#'   \code{get.spChkOption()}.
#' @param adjust.n Default TRUE, if FALSE the number of
#'   observations is not adjusted for no-neighbor
#'   observations, if TRUE, the number of observations is
#'   adjusted.
#' @param B1correct Default TRUE, if TRUE, the erratum
#'   referenced below: "On page 195, the coefficient of W2
#'   in B1, (just below center of the page) should be 6,
#'   not 3." is applied; if FALSE, 3 is used (as in CrimeStat
#'   IV).
#' @param adjust.x Default TRUE, if TRUE, x values of
#'   observations with no neighbors are omitted in the
#'   denominator of G.
#' @param Arc_all_x Default FALSE, if Arc_all_x=TRUE and
#'   adjust.x=TRUE, use the full x vector in part of the
#'   denominator term for G.
#' @return A list with class \code{htest} containing the
#'   following components:
#'   \item{statistic}{the value of the standard deviate
#'     of Getis-Ord G.}
#'   \item{p.value}{the p-value of the test.}
#'   \item{estimate}{the value of the observed statistic,
#'     its expectation and variance.}
#'   \item{alternative}{a character string describing the
#'     alternative hypothesis.}
#'   \item{data.name}{a character string giving the
#'     name(s) of the data.}
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
getisGlobalGTest <- function(x,
                             listw,
                             zero.policy = NULL,
                             alternative = "greater",
                             spChk = NULL,
                             adjust.n = TRUE,
                             B1correct = TRUE,
                             adjust.x = TRUE,
                             Arc_all_x = FALSE) {
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
#' @param x Numeric vector the same length as the neighbors
#' list in listw
#' @param listw A \code{listw} object created, for example,
#' by \code{nb2listw}
#' @param zero.policy Default NULL, use global option value;
#' if TRUE assign zero to the lagged value of zones without
#' neighbors, if FALSE assign NA
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
#' @return A vector of G or Gstar standard deviate values,
#' with attributes "gstari" set to TRUE or FALSE, "call"
#' set to the function call, and class "localG". For
#' conditional permutation, the returned value is the same
#' as for \code{localG()}, and the simulated standard deviate
#' is returned as column \code{"StdDev.Gi"} in
#' \code{attr(., "internals")}.
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
getisLocalG <- function(x,
                       listw,
                       zero.policy = NULL,
                       spChk = NULL,
                       GeoDa = FALSE,
                       alternative = "two.sided") {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localG function from spdep
  spdep::localG(x = x,
                listw = listw,
                zero.policy = zero.policy,
                spChk = spChk,
                GeoDa = GeoDa,
                alternative = alternative,
                return_internals = return_internals)
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
getisLocalGPerm <- function(x,
                            listw,
                            nsim = 499,
                            zero.policy = NULL,
                            spChk = NULL,
                            alternative = "two.sided",
                            iseed = NULL,
                            fix_i_in_Gstar_permutations = TRUE,
                            no_repeat_in_row = FALSE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localG_perm function from spdep
  spdep::localG_perm(x = x,
                     listw = listw,
                     nsim = nsim,
                     zero.policy = zero.policy,
                     spChk = spChk,
                     alternative = alternative,
                     iseed = iseed,
                     fix_i_in_Gstar_permutations = fix_i_in_Gstar_permutations,
                     no_repeat_in_row = no_repeat_in_row)
}
