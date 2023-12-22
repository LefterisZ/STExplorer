#' Compute Moran's I
#'
#' This function is a wrapper of the \code{\link[spdep]{moran}} function from
#' the `spdep` package developed by Roger Bivand. It is a function to compute
#' Moran's I, called by \code{moran.test} and \code{moran.mc}.
#'
#' @param x a numeric vector same length as neighbours in listw
#' @param listw \code{listw} object, e.g., by \code{nb2listw}
#' @param n number of zones
#' @param S0 global sum of weights
#' @param zero.policy default \code{attr(listw, "zero.policy")}
#' @param NAOK if 'TRUE' pass 'NA' or 'NaN' or 'Inf' to foreign function
#'
#' @return list of
#'  \item{I}{Moran's I}
#'  \item{K}{sample kurtosis of x}
#'
#' @references Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 17.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname moranGlobalI
#'
#' @aliases moranGlobalI
#'
#' @seealso \code{\link{moran.test}}, \code{\link{moran.mc}}
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
#' moran <- moranGlobalI(x = gene_exp, listw = listw,zero.policy = TRUE)
#' }
#'
#' @export
moranGlobalI <- function(x,
                         listw,
                         n = NULL,
                         S0 = NULL,
                         zero.policy = NULL,
                         NAOK = FALSE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Calculate n and S0 if not provided
  if (is.null(n)) n <- length(listw$neighbours)
  if (is.null(S0)) S0 <- spdep::Szero(listw)

  ## Call the moran function from spdep
  spdep::moran(x,
               listw,
               n = n,
               S0 = S0,
               zero.policy = zero.policy,
               NAOK = NAOK)
}


#' Permutation test for Moran's I statistic
#'
#' This function is a wrapper of the \code{\link[spdep]{moran.mc}} function from
#' the `spdep` package developed by Roger Bivand. A permutation test for
#' Moran's I statistic calculated by using nsim random
#' permutations of x for the given spatial weighting scheme, establishing the
#' rank of the observed statistic in relation to the nsim simulated values.
#'
#' @param x Numeric vector, same length as neighbours list in listw
#' @param listw \code{listw} object, e.g., from \code{nb2listw}
#' @param nsim Number of permutations
#' @param zero.policy Default \code{attr(listw, "zero.policy")} as set in
#' \code{listw}
#' @param alternative Character string specifying the alternative hypothesis.
#' Options: "greater" (default), "two.sided", or "less".
#' @param na.action Function (default \code{na.fail}). Options: \code{na.omit},
#' \code{na.exclude}. The weights list is subsetted to remove NAs in the data.
#' \code{na.pass} is not permitted because it is meaningless in a permutation
#' test.
#' @param spChk Check data vector names against spatial objects for identity
#' integrity. TRUE or FALSE, default NULL to use \code{get.spChkOption()}
#' @param return_boot Return an object of class \code{boot} from the equivalent
#' permutation bootstrap.
#' @param adjust.n Default TRUE. If FALSE, the number of observations is not
#' adjusted for no-neighbour observations. If TRUE, the number of observations
#' is adjusted.
#'
#' @return A list with class \code{htest} and \code{mc.sim} containing the
#' following components:
#' \item{statistic}{the value of the observed Moran's I.}
#' \item{parameter}{the rank of the observed Moran's I.}
#' \item{p.value}{the pseudo p-value of the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string giving the method used.}
#' \item{data.name}{a character string giving the name(s) of the data, and the
#' number of simulations.}
#' \item{res}{nsim simulated values of the statistic, final value is observed
#' statistic}
#'
#' @references Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.
#'
#' @seealso \code{\link{moran}}, \code{\link{moran.test}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname moranGlobalIPerm
#'
#' @aliases moranGlobalIPerm
#'
#' @seealso \code{\link{moran.test}}, \code{\link{moran}}
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
#' moran <- moranGlobalIPerm(x = gene_exp, listw = listw, nsim = 499,
#' zero.policy = TRUE)
#' }
#'
#' @export
moranGlobalIPerm <- function(x,
                             listw,
                             nsim,
                             zero.policy = NULL,
                             alternative = "greater",
                             na.action = na.fail,
                             spChk = NULL,
                             return_boot = FALSE,
                             adjust.n = TRUE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the moran.mc function from spdep
  spdep::moran.mc(x = x,
                  listw = listw,
                  nsim = nsim,
                  zero.policy = zero.policy,
                  alternative = alternative,
                  na.action = na.action,
                  spChk = spChk,
                  return_boot = return_boot,
                  adjust.n = adjust.n)
}


#' Moran's I test for spatial autocorrelation
#'
#' This function is a wrapper of the \code{\link[spdep]{moran.mc}} function from
#' the `spdep` package developed by Roger Bivand. Moran's test for spatial
#' autocorrelation using a spatial weights matrix in weights list form. The
#' assumptions underlying the test are sensitive to the form of the graph of
#' neighbour relationships and other factors, and results may be checked
#' against those of \code{moran.mc} permutations.
#'
#' @param x Numeric vector, same length as neighbors list in listw
#' @param listw \code{listw} object, e.g., from \code{nb2listw}
#' @param randomisation Variance of I calculated under the assumption of
#' randomization. If FALSE, normality.
#' @param zero.policy Default \code{attr(listw, "zero.policy")} as set in
#' \code{listw}
#' @param alternative Character string specifying the alternative hypothesis.
#' Options: greater (default), less, or two.sided.
#' @param rank Logical value - default FALSE for continuous variables.
#' If TRUE, uses the adaptation of Moran's I for ranks suggested by Cliff and
#' Ord (1981, p. 46).
#' @param na.action Function (default \code{na.fail}). Options: \code{na.omit},
#' \code{na.exclude}. The weights list is subsetted to remove NAs in the data.
#' \code{na.pass} is used, zero is substituted for NA values in calculating the
#' spatial lag.
#' @param spChk Check data vector names against spatial objects for identity
#' integrity.
#' TRUE or FALSE, default NULL to use \code{get.spChkOption()}
#' @param adjust.n Default TRUE. If FALSE, the number of observations is not
#' adjusted for no-neighbor observations. If TRUE, the number of observations
#' is adjusted.
#' @param drop.EI2 Default FALSE. If TRUE, emulate CrimeStat <= 4.02.
#'
#' @return A list with class \code{htest} containing the following components:
#' \item{statistic}{the value of the standard deviate of Moran's I.}
#' \item{p.value}{the p-value of the test.}
#' \item{estimate}{the value of the observed Moran's I, its expectation and
#' variance under the method assumption.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string giving the assumption used for calculating
#' the standard deviate.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#'
#' @note Var(I) is taken from Cliff and Ord (1969, p. 28),
#' and Goodchild's CATMOG 47 (1986), see also Upton & Fingleton (1985) p. 171;
#' it agrees with SpaceStat, see Tutorial workbook Chapter 22; VI is the second
#' crude moment minus the square of the first crude moment. The derivation of the
#' test (Cliff and Ord, 1981, p. 18) assumes that the weights matrix is symmetric.
#' For inherently non-symmetric matrices, such as k-nearest neighbor matrices,
#' \code{listw2U()} can be used to make the matrix symmetric.
#'
#' @references Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21;
#' Bivand RS, Wong DWS 2018 Comparing implementations of global and local
#' indicators of spatial association. TEST, 27(3), 716--748 \doi{10.1007/s11749-
#' 018-0599-x}
#'
#' @seealso \code{\link{moran}}, \code{\link{moran.mc}}, \code{\link{listw2U}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname moranGlobalITest
#'
#' @aliases moranGlobalITest
#'
#' @seealso \code{\link{moran}}, \code{\link{moran.mc}}
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
#' # A listw object containing distance-based spatial weights for neighbors
#' listw <- listw
#'
#' # Calculate Moran's I for this gene
#' moran <- moranGlobalITest(x = gene_exp, listw = listw, zero.policy = TRUE)
#' }
#'
#' @export
moranGlobalITest <- function(x,
                             listw,
                             randomisation = TRUE,
                             zero.policy = NULL,
                             alternative = "greater",
                             rank = FALSE,
                             na.action = na.fail,
                             spChk = NULL,
                             adjust.n = TRUE,
                             drop.EI2 = FALSE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the moran.test function from spdep
  spdep::moran.test(x = x,
                    listw = listw,
                    randomisation = randomisation,
                    zero.policy = zero.policy,
                    alternative = alternative,
                    rank = rank,
                    na.action = na.action,
                    spChk = spChk,
                    adjust.n = adjust.n,
                    drop.EI2 = drop.EI2)
}


#' Local Moran's I statistic
#'
#' This function is a wrapper of the \code{\link[spdep]{moran}} function from
#' the `spdep` package developed by Roger Bivand. The local spatial statistic
#' Moran's I is calculated for each zone based on the spatial weights object
#' used. The values returned include a Z-value, and may be used as a diagnostic
#' tool.
#'
#' @param x Numeric vector the same length as the neighbors
#' list in listw
#' @param listw \code{listw} object created, for example,
#' by \code{nb2listw}
#' @param zero.policy Default \code{attr(listw,
#' "zero.policy")} as set when \code{listw} was created. If
#' attribute not set, use global option value. If TRUE,
#' assign zero to the lagged value of zones without
#' neighbors; if FALSE, assign NA.
#' @param na.action A function (default \code{na.fail}), can
#' also be \code{na.omit} or \code{na.exclude} - in these
#' cases the weights list will be subsetted to remove NAs
#' in the data. It may be necessary to set zero.policy to
#' TRUE because this subsetting may create no-neighbor
#' observations. Note that only weights lists created
#' without using the glist argument to \code{nb2listw} may
#' be subsetted. If \code{na.pass} is used, zero is
#' substituted for NA values in calculating the spatial
#' lag. (Note that na.exclude will only work properly
#' starting from R 1.9.0, na.omit and na.exclude assign the
#' wrong classes in 1.8.*)
#' @param conditional Default TRUE: expectation and variance
#' are calculated using the conditional randomization null
#' (Sokal 1998 Eqs. A7 & A8). Elaboration of these changes
#' available in Sauer et al. (2021). If FALSE: expectation
#' and variance are calculated using the total
#' randomization null (Sokal 1998 Eqs. A3 & A4).
#' @param alternative A character string specifying the
#' alternative hypothesis, must be one of greater, less or
#' two.sided (default).
#' @param mlvar Default TRUE: values of local Moran's I are
#' reported using the variance of the variable of interest
#' (sum of squared deviances over n), but can be reported
#' as the sample variance, dividing by (n-1) instead; both
#' are used in other implementations.
#' @param spChk Should the data vector names be checked
#' against the spatial objects for identity integrity,
#' TRUE, or FALSE, default NULL to use
#' \code{get.spChkOption()}.
#' @param adjust.x Default FALSE, if TRUE, x values of
#' observations with no neighbors are omitted in the mean
#' of x.
#' @param nsim Default 499, number of conditional
#' permutation simulations.
#' @param sample_Ei Default TRUE; if conditional permutation,
#' use the sample $E_i$ values, or the analytical values,
#' leaving only variances calculated by simulation.
#' @param iseed Default NULL, used to set the seed for
#' possible parallel RNGs.
#' @param no_repeat_in_row Default \code{FALSE}, if
#' \code{TRUE}, sample conditionally in each row without
#' replacements to avoid duplicate values,
#' \url{https://github.com/r-spatial/spdep/issues/124}.

#'
#' @details The values of local Moran's I are divided by the variance (or
#' sample variance) of the variable of interest to accord with Table 1, p. 103,
#' and formula (12), p. 99, in Anselin (1995), rather than his formula (7), p.
#' 98. The variance of the local Moran statistic is taken from Sokal et al.
#' (1998) p. 334, equations 4 & 5 or equations 7 & 8 located depending on user
#' specification. By default, the implementation divides by n, not (n-1) in
#' calculating the variance and higher moments. Conditional code contributed by
#' Jeff Sauer and Levi Wolf.
#'
#' @return A list with the following components:
#' \item{Ii}{local Moran statistic}
#' \item{E.Ii}{expectation of local Moran statistic; for \code{localmoran_perm}
#' the permutation sample means}
#' \item{Var.Ii}{variance of local Moran statistic; for \code{localmoran_perm}
#' the permutation sample standard deviations}
#' \item{Z.Ii}{standard deviate of local Moran statistic; for
#' \code{localmoran_perm} based on permutation sample means and standard
#' deviations}
#' \item{Pr()}{p-value of local Moran statistic using \code{pnorm()}; for
#' \code{localmoran_perm} using standard deviates based on permutation sample
#' means and standard deviations}
#' \item{Pr() Sim}{For \code{localmoran_perm}, \code{rank()} and \code{punif()}
#' of observed statistic rank for [0, 1] p-values using \code{alternative=}}
#' \item{Pr(folded) Sim}{the simulation folded [0, 0.5] range ranked p-value
#' (based on \url{https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b
#' 8cadcbecc5e061/esda/crand.py#L211-L213})}
#' \item{Skewness}{For \code{localmoran_perm}, the output of
#' \code{e1071::skewness()} for the permutation samples underlying the standard
#' deviates}
#' \item{Kurtosis}{For \code{localmoran_perm}, the output of
#' \code{e1071::kurtosis()} for the permutation samples underlying the standard
#' deviates}
#'
#' In addition, an attribute data frame \code{"quadr"} with mean and median
#' quadrant columns, and a column splitting on the demeaned variable and lagged
#' demeaned variable at zero.
#'
#' @note Conditional permutations added for comparative purposes; permutations
#' are over the whole data vector omitting the observation itself. For p-value
#' adjustment, use \code{p.adjust()} or \code{p.adjustSP()} on the output
#' vector.
#'
#' @references Anselin, L. 1995. Local indicators of spatial association,
#' Geographical Analysis, 27, 93--115;
#' Getis, A. and Ord, J. K. 1996 Local spatial statistics: an overview. In P.
#' Longley and M. Batty (eds) \emph{Spatial analysis: modelling in a GIS
#' environment} (Cambridge: Geoinformation International), 261--277;
#' Sokal, R. R, Oden, N. L. and Thomson, B. A. 1998. Local Spatial
#' Autocorrelation in a Biological Model. Geographical Analysis, 30. 331--354;
#' Bivand RS, Wong DWS 2018 Comparing implementations of global and local
#' indicators of spatial association. TEST, 27(3), 716--748 \doi{10.1007/s11749-
#' 018-0599-x};
#' Sauer, J., Oshan, T. M., Rey, S., & Wolf, L. J. 2021. The Importance of Null
#' Hypotheses: Understanding Differences in Local Moran’s under
#' Heteroskedasticity. Geographical Analysis. \doi{doi:10.1111/gean.12304}
#' Bivand, R. (2022), R Packages for Analyzing Spatial Data: A Comparative Case
#' Study with Areal Data. Geographical Analysis, 54(3), 488-518. \doi{10.1111/
#' gean.12319}
#'
#' @seealso \code{\link{localG}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname moranLocalI
#'
#' @aliases moranLocalI
#'
#' @seealso \code{\link{localmoran_perm}}
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
#' # Calculate local Moran's I for this gene
#' moran <- moranLocalI(x = gene_exp, listw = listw,zero.policy = TRUE)
#' }
#'
#' @export
moranLocalI <- function(x,
                        listw,
                        zero.policy = NULL,
                        na.action = na.fail,
                        conditional = TRUE,
                        alternative = "two.sided",
                        mlvar = TRUE,
                        spChk = NULL,
                        adjust.x = FALSE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localmoran function from spdep
  spdep::localmoran(x = x,
                    listw = listw,
                    zero.policy = zero.policy,
                    na.action = na.action,
                    conditional = conditional,
                    alternative = alternative,
                    mlvar = mlvar,
                    spChk = spChk,
                    adjust.x = adjust.x)
}


#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname moranLocalI
#'
#' @aliases moranLocalIPerm
#'
#' @examples
#' \dontrun{
#' # Calculate local Moran's I for this gene
#' moran <- moranLocalIPerm(x = gene_exp, listw = listw, nsim = 499,
#' zero.policy = TRUE)
#' }
#'
#' @export
moranLocalIPerm <- function(x,
                            listw,
                            nsim = 499,
                            zero.policy = NULL,
                            na.action = na.fail,
                            alternative = "two.sided",
                            mlvar = TRUE,
                            spChk = NULL,
                            adjust.x = FALSE,
                            sample_Ei = TRUE,
                            iseed = NULL,
                            no_repeat_in_row = FALSE) {
  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localmoran_perm function from spdep
  spdep::localmoran_perm(x = x,
                         listw = listw,
                         nsim = nsim,
                         zero.policy = zero.policy,
                         na.action = na.action,
                         alternative = alternative,
                         mlvar = mlvar,
                         spChk = spChk,
                         adjust.x = adjust.x,
                         sample_Ei = sample_Ei,
                         iseed = iseed,
                         no_repeat_in_row = no_repeat_in_row)
}


