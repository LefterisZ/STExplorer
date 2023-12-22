#' Calculate Geary's C Global Statistic
#'
#' This function is a wrapper of the \code{\link[spdep]{geary}} function from
#' the `spdep` package developed by Roger Bivand. It computes Geary's C
#' statistic for spatial autocorrelation using the spdep package.
#'
#' @param x A numeric vector of observations. Must be the same length as the
#' neighbours list in listw.
#' @param listw A listw object created, for example, by
#' \code{link[spdep]{nb2listw}} or \code{link[spdep]{nb2listwdist}}.
#' @param n Number of zones. If not provided, it is calculated as
#' length(listw$neighbours).
#' @param n1 \code{n - 1}. If not provided, it is calculated as
#' \code{length(listw$neighbours) - 1}.
#' @param S0 Global sum of weights. If not provided, it is calculated as
#' \code{spdep::Szero(listw)}.
#' @param zero.policy Default \code{attr(listw, "zero.policy")} as set when
#' \code{listw} was created.
#'   If TRUE, assign zero to the lagged value of zones without neighbours.
#'   If FALSE, assign NA. If attribute not set, use the global option value.
#'
#' @return Geary's C statistic. A list with C: Geary's C and K: sample kurtosis
#' of x
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
gearyGlobalC <- function(x,
                         listw,
                         n = NULL,
                         n1 = NULL,
                         S0 = NULL,
                         zero.policy = NULL) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Calculate n, n1, and S0 if not provided
  if (is.null(n)) n <- length(listw$neighbours)
  if (is.null(n1)) n1 <- n - 1
  if (is.null(S0)) S0 <- spdep::Szero(listw)

  ## Call the geary function from spdep
  spdep::geary(x, listw, n = n, n1 = n1, S0 = S0, zero.policy = zero.policy)
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
#' @param x A numeric vector of observations. Must be the same length as the
#' neighbours list in listw.
#' @param listw A listw object created, for example, by
#' \code{link[spdep]{nb2listw}} or \code{link[spdep]{nb2listwdist}}.
#' @param nsim Number of permutations.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "greater" (default), or "less"; this reversal corresponds to
#' that on geary.test described in the section on the output statistic value,
#' based on Cliff and Ord 1973, p. 21.
#' @param spChk should the data vector names be checked against the spatial
#' objects for identity integrity, TRUE, or FALSE, default NULL to use
#' \code{get.spChkOption()}.
#' @param zero.policy Default \code{attr(listw, "zero.policy")} as set when
#' \code{listw} was created.
#'   If TRUE, assign zero to the lagged value of zones without neighbours.
#'   If FALSE, assign NA. If attribute not set, use the global option value.
#' @param adjust.n Default TRUE, if FALSE the number of observations is not
#' adjusted for no-neighbour observations, if TRUE, the number of observations
#' is adjusted.
#' @param return_boot Return an object of class boot from the equivalent
#' permutation bootstrap rather than an object of class htest
#'
#'
#' @return A list with class \code{htest} and \code{mc.sim} containing the
#' following components:
#' \itemize{
#'  \item{statistic}{the value of the observed Geary's C.}
#'  \item{parameter}{the rank of the observed Geary's C.}
#'  \item{p.value}{the pseudo p-value of the test.}
#'  \item{alternative}{a character string describing the alternative hypothesis.}
#'  \item{method}{a character string giving the method used.}
#'  \item{data.name}{a character string giving the name(s) of the data, and the
#'  number of simulations.}
#'  \item{res}{nsim simulated values of statistic, final value is observed
#'  statistic}
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
gearyGlobalCPerm <- function(x,
                             listw,
                             nsim,
                             zero.policy = attr(listw, "zero.policy"),
                             alternative = "greater",
                             spChk = NULL,
                             adjust.n = TRUE,
                             return_boot = FALSE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Assign nsim if not provided
  if (missing(nsim)) nsim <- 999

  ## Call the geary.mc function from spdep
  spdep::geary.mc(x = x,
                  listw = listw,
                  nsim = nsim,
                  alternative = alternative,
                  spChk = spChk,
                  adjust.n = adjust.n,
                  return_boot = return_boot,
                  zero.policy = zero.policy)
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
#' @param x A numeric vector the same length as the neighbours list in listw.
#' @param listw A \code{listw} object created, for example, by \code{nb2listw}.
#' @param randomisation Variance of I calculated under the assumption of
#' randomisation. If FALSE, normality.
#' @param zero.policy Default \code{attr(listw, "zero.policy")} as set when
#' \code{listw} was created. If attribute not set, use the global option value.
#' If TRUE, assign zero to the lagged value of zones without neighbours, if
#' FALSE, assign NA.
#' @param alternative A character string specifying the alternative hypothesis.
#' Must be one of "greater" (default), "less" or "two.sided".
#' @param spChk Should the data vector names be checked against the spatial
#' objects for identity integrity? TRUE or FALSE, default NULL to use
#' \code{get.spChkOption()}.
#' @param adjust.n Default TRUE. If FALSE, the number of observations is not
#' adjusted for no-neighbour observations. If TRUE, the number of observations
#' is adjusted.
#'
#' @return A list with class \code{htest} containing the following components:
#'   \item{statistic}{the value of the standard deviate of Geary's C.}
#'   \item{p.value}{the p-value of the test.}
#'   \item{estimate}{the value of the observed Geary's C, its expectation, and
#'   variance under the method assumption.}
#'   \item{alternative}{a character string describing the alternative
#'   hypothesis.}
#'   \item{method}{a character string giving the assumption used for
#'   calculating the standard deviate.}
#'   \item{data.name}{a character string giving the name(s) of the data.}
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
gearyGlobalCTest <- function(x,
                             listw,
                             randomisation=TRUE,
                             zero.policy = attr(listw, "zero.policy"),
                             alternative = "greater",
                             spChk = NULL,
                             adjust.n = TRUE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the geary.test function from spdep
  spdep::geary.test(x = x,
                    listw = listw,
                    randomisation = randomisation,
                    alternative = alternative,
                    spChk = spChk,
                    adjust.n = adjust.n,
                    zero.policy = zero.policy)
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
#' @param x A numeric vector, list, or matrix. In univariate
#' context, \code{x} can be a numeric vector. In multivariate
#' context, provide a list or matrix.
#' @param listw A \code{listw} object created, e.g., by
#' \code{nb2listw}.
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
#' @returns Numeric vector containing Local Geary statistic
#' with attribute \code{pseudo-p} when \code{localC_perm()} is
#' used. \code{pseudo-p} is an 8-column matrix containing:
#'   \item{E.Ci}{expectation of Local Geary statistic based on
#' permutation sample}
#'   \item{Var.Ci}{variance of Local Geary based on permutation
#' sample}
#'   \item{Z.Ci}{standard deviate of Local Geary based on
#' permutation sample}
#'   \item{Pr()}{p-value of Local Geary statistic using
#' \code{pnorm()} based on permutation sample means and standard
#' deviations}
#'   \item{Pr() Sim}{\code{rank()} and \code{punif()} of
#' observed statistic rank for [0, 1] p-values using
#' \code{alternative=}}
#'   \item{Pr(folded) Sim}{simulation folded [0, 0.5] range
#' ranked p-value (based on
#' \url{https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b8cadcbecc
#' 5e061/esda/crand.py#L211-L213})}
#'   \item{Skewness}{output of \code{e1071::skewness()} for
#' permutation samples underlying standard deviates}
#'   \item{Kurtosis}{output of \code{e1071::kurtosis()} for
#' permutation samples underlying standard deviates}
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
gearyLocalC <- function(x,
                        listw,
                        ...,
                        zero.policy = attr(listw, "zero.policy")) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  spdep::localC(x = x,
                listw = listw,
                ...,
                zero.policy = zero.policy)
}


#' @aliases gearyLocalCPerm
#' @rdname gearyLocalC
#' @examples
#' \dontrun{
#' # Calculate  Geary's C and statistical significance for this gene
#' geary <- gearyLocalCPerm(x = gene_exp, listw = listw, nsim = 499,
#' zero.policy = TRUE)
#' }
#' @export
gearyLocalCPerm <- function(x,
                            listw,
                            nsim = 499,
                            alternative = "two.sided",
                            ...,
                            zero.policy = attr(listw, "zero.policy"),
                            iseed = NULL,
                            no_repeat_in_row = FALSE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  spdep::localC_perm(x = x,
                     listw = listw,
                     nsim = nsim,
                     alternative = alternative,
                     ...,
                     zero.policy = zero.policy,
                     iseed = iseed,
                     no_repeat_in_row = no_repeat_in_row)

}


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
#'
#' @export
.int_checkSAInput <- function(x, listw) {
  if (!is.numeric(x) ||
      length(x) != length(listw$neighbours) ||
      !inherits(listw, "listw")) {
    stop("Invalid input. Please provide a numeric vector 'x', listw",
         " of class 'listw', and ensure the length of 'x' matches the",
         " length of listw$neighbours.")
  }
}
