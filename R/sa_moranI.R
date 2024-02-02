#' Compute Moran's I
#'
#' This function is a wrapper of the \code{\link[spdep]{moran}} function from
#' the `spdep` package developed by Roger Bivand. It is a function to compute
#' Moran's I, called by \code{moran.test} and \code{moran.mc}.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
#' @param n number of zones
#' @param S0 global sum of weights
#' @param NAOK if 'TRUE' pass 'NA' or 'NaN' or 'Inf' to foreign function
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
#'  \item{moranI}{Moran's I}
#'  \item{moranK}{sample kurtosis of x}
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
moranGlobalI <- function(m_sfe,
                         sample_id,
                         genes = TRUE,
                         n = NULL,
                         S0 = NULL,
                         zero.policy = NULL,
                         NAOK = FALSE,
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

  ## Calculate n and S0 if not provided
  if (is.null(n)) n <- length(listw$neighbours)
  if (is.null(S0)) S0 <- spdep::Szero(listw)

  ## Calculate global moran's I
  out <- parallel::mclapply(genes,
                            .int_moran,
                            sfe = sfe,
                            listw = listw,
                            n = n,
                            S0 = S0,
                            NAOK = NAOK,
                            zero.policy = zero.policy,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's rowData
  out <- as.data.frame(rlist::list.rbind(out))
  SummarizedExperiment::rowData(sfe)$moranI <- out$I
  SummarizedExperiment::rowData(sfe)$moranK <- out$K

  return(sfe)
}


#' Permutation test for Moran's I statistic
#'
#' This function is a wrapper of the \code{\link[spdep]{moran.mc}} function from
#' the `spdep` package developed by Roger Bivand. A permutation test for
#' Moran's I statistic calculated by using nsim random
#' permutations of x for the given spatial weighting scheme, establishing the
#' rank of the observed statistic in relation to the nsim simulated values.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
#' @param nsim Number of permutations
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
#' @return An SFE object with the below two columns added in the `rowData`:
#' \item{I}{the value of the observed Moran's I. - replaces other global I
#' calculations since the I is always going to be the same.}
#' \item{moranPval_perm}{the pseudo p-value of the test.}
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
moranGlobalIPerm <- function(m_sfe,
                             sample_id,
                             genes = TRUE,
                             nsim = 999,
                             zero.policy = NULL,
                             alternative = "greater",
                             na.action = na.fail,
                             spChk = NULL,
                             return_boot = FALSE,
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

  out <- parallel::mclapply(genes,
                            .int_moranIPerm,
                            sfe = sfe,
                            listw = listw,
                            nsim = nsim,
                            alternative = alternative,
                            na.action = na.action,
                            spChk = spChk,
                            return_boot = return_boot,
                            adjust.n = adjust.n,
                            zero.policy = zero.policy,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's rowData
  out <- as.data.frame(rlist::list.rbind(out))
  SummarizedExperiment::rowData(sfe)$moranI <- unlist(out$statistic)
  SummarizedExperiment::rowData(sfe)$moranPval_perm <- unlist(out$p.value)

  return(sfe)
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
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#' @param sample_id A character vector specifying the sample IDs to include
#' (only relevant if a MetaSpatialFeatureExperiment has been provided in the
#' `m_sfe` argument).
#' @param genes TRUE or a named character vector with gene names for which the
#' SA statistic needs to be calculated. If left to TRUE then the SA statistic
#' is calculated for every gene.
#' @param randomisation Variance of I calculated under the assumption of
#' randomization. If FALSE, normality.
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
#' @return An SFE object with the below two columns added in the `rowData`:
#' \item{moranI_test}{the value of the observed Moran's I.}
#' \item{moranPval_test}{the pseudo p-value of the test.}
#'
#' @note Var(I) is taken from Cliff and Ord (1969, p. 28),
#' and Goodchild's CATMOG 47 (1986), see also Upton & Fingleton (1985) p. 171;
#' it agrees with SpaceStat, see Tutorial workbook Chapter 22; VI is the second
#' crude moment minus the square of the first crude moment. The derivation of
#' the test (Cliff and Ord, 1981, p. 18) assumes that the weights matrix is
#' symmetric. For inherently non-symmetric matrices, such as k-nearest
#' neighbour matrices, \code{listw2U()} can be used to make the matrix
#' symmetric.
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
moranGlobalITest <- function(m_sfe,
                             sample_id,
                             genes = TRUE,
                             zero.policy = NULL,
                             randomisation = TRUE,
                             alternative = "greater",
                             rank = FALSE,
                             na.action = na.fail,
                             spChk = NULL,
                             adjust.n = TRUE,
                             drop.EI2 = FALSE,
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
                            .int_moranITest,
                            sfe = sfe,
                            listw = listw,
                            randomisation = randomisation,
                            alternative = alternative,
                            rank = rank,
                            na.action = na.action,
                            spChk = spChk,
                            adjust.n = adjust.n,
                            drop.EI2 = drop.EI2,
                            zero.policy = zero.policy,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's rowData
  out <- as.data.frame(rlist::list.rbind(out))
  SummarizedExperiment::rowData(sfe)$moranI <- unlist(out$statistic)
  SummarizedExperiment::rowData(sfe)$moranPval_test <- unlist(out$p.value)

  return(sfe)
}


#' Local Moran's I statistic
#'
#' This function is a wrapper of the \code{\link[spdep]{moran}} function from
#' the `spdep` package developed by Roger Bivand. The local spatial statistic
#' Moran's I is calculated for each zone based on the spatial weights object
#' used. The values returned include a Z-value, and may be used as a diagnostic
#' tool.
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
#' @return An SFE object with the results added in the `localResults` slot
#' which contains a DataFrame named `localMoranI`/`localMoranIPerm` that
#' contains the per location Ii statistics values for each gene:
#' \item{ENSG***.Ii}{Numeric vector containing Local Moran statistic
#' with attribute
#' \item{ENSG***.IiFDR}{(Only for gearyLocalCPerm) \code{rank()} and
#' \code{punif()} of observed statistic rank for [0, 1] p-values using
#' \code{alternative=}}
#' \item{ENSG***.IiCluster}{Low/High cluster for each location. It can take
#' values: "Low-Low", "High-High", "Low-High", "High-Low"}
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
moranLocalI <- function(m_sfe,
                        sample_id,
                        genes = TRUE,
                        zero.policy = NULL,
                        na.action = na.fail,
                        conditional = TRUE,
                        alternative = "two.sided",
                        mlvar = TRUE,
                        spChk = NULL,
                        adjust.x = FALSE,
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
                            .int_moranLocalI,
                            sfe = sfe,
                            listw = listw,
                            na.action = na.action,
                            conditional = conditional,
                            alternative = alternative,
                            mlvar = mlvar,
                            spChk = spChk,
                            adjust.x = adjust.x,
                            zero.policy = zero.policy,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's localResults
  out <- S4Vectors::DataFrame(rlist::list.cbind(out))
  localResults(sfe, name = "localMoranI") <- out

  return(sfe)
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
moranLocalIPerm <- function(m_sfe,
                            sample_id,
                            genes = TRUE,
                            nsim = 999,
                            zero.policy = NULL,
                            na.action = na.fail,
                            alternative = "two.sided",
                            mlvar = TRUE,
                            spChk = NULL,
                            adjust.x = FALSE,
                            sample_Ei = TRUE,
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
    stop("Invalid `genes` argument input")
  }

  ## Get neighbour graph
  listw <- colGraph(sfe)

  ## Check zero.policy
  if (is.null(zero.policy)) {
    zero.policy = attr(listw, "zero.policy")
  }

  out <- parallel::mclapply(genes,
                            .int_moranLocalIPerm,
                            sfe = sfe,
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
                            no_repeat_in_row = no_repeat_in_row,
                            mc.cores = mc.cores)

  ## Import output into the SFE object's localResults
  out <- S4Vectors::DataFrame(rlist::list.cbind(out))
  localResults(sfe, name = "localMoranIPerm") <- out

  return(sfe)
}


# ---------------------------------------------------------------------------- #
#  ###### INTERNAL FUNCTIONS ASSOCIATED WITH SA CALCULATIONS (C, G, I) ######
# ---------------------------------------------------------------------------- #
#' Internal Function: Calculate Moran's I
#'
#' This internal function calculates the I statistic for a vector of gene
#' expression. For the parameter arguments check the
#' \code{link[STExplorer]{moranGlobalI}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @importFrom SummarizedExperiment assay
#'
#' @rdname dot-int_moran
#'
#' @aliases .int_moran
.int_moran <- function(gene,
                       sfe,
                       listw,
                       n,
                       S0,
                       zero.policy,
                       NAOK) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the moran function from spdep
  spdep::moran(x,
               listw,
               n = n,
               S0 = S0,
               zero.policy = zero.policy,
               NAOK = NAOK)
}


#' Internal Function: Calculate Moran's I with Permutation testing
#'
#' This internal function calculates the I statistic for a vector of gene
#' expression with Monte-Carlo permutation testing of significance. For the
#' parameter arguments check the \code{link[STExplorer]{moranGlobalIPerm}}
#' function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_moranIPerm
#'
#' @aliases .int_moranIPerm
#'
.int_moranIPerm <- function(gene,
                            sfe,
                            listw,
                            nsim,
                            alternative,
                            na.action,
                            spChk,
                            adjust.n,
                            return_boot,
                            zero.policy) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

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


#' Internal Function: Calculate Moran's I with statistical testing
#'
#' This internal function calculates the I statistic for a vector of gene
#' expression alongside a statistical significance testing. For the parameter
#' arguments check the \code{link[STExplorer]{moranGlobalITest}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_moranITest
#'
#' @aliases .int_moranITest
#'
.int_moranITest <- function(gene,
                            sfe,
                            listw,
                            randomisation,
                            alternative,
                            rank,
                            na.action,
                            spChk,
                            adjust.n,
                            drop.EI2,
                            zero.policy) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

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


#' Internal Function: Calculate local Moran's I
#'
#' This internal function calculates the local Ii statistic for a vector of
#' gene expression. For the parameter arguments check the
#' \code{link[STExplorer]{moranLocalI}} function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_moranLocalI
#'
#' @aliases .int_moranLocalI
#'
.int_moranLocalI <- function(gene,
                             sfe,
                             listw,
                             zero.policy,
                             na.action,
                             conditional,
                             alternative,
                             mlvar,
                             spChk,
                             adjust.x) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localmoran function from spdep
  res <- spdep::localmoran(x = x,
                           listw = listw,
                           zero.policy = zero.policy,
                           na.action = na.action,
                           conditional = conditional,
                           alternative = alternative,
                           mlvar = mlvar,
                           spChk = spChk,
                           adjust.x = adjust.x)

  ## Get attributes into a dataframe
  out <- data.frame(Ii = res[,"Ii"],
                    IiFDR = res[,dimnames(res)[[2]][5]],
                    IiClust = attr(res, "quadr")[,"pysal"]) # maybe give option to select between mean, median, and pysal

  return(out)
}


#' Internal Function: Calculate local Moran's I with statistical testing
#'
#' This internal function calculates the local Ii statistic for a vector of
#' gene expression alongside a statistical significance testing. For the
#' parameter arguments check the \code{link[STExplorer]{moranLocalIPerm}}
#' function.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_moranLocalIPerm
#'
#' @aliases .int_moranLocalIPerm
#'
.int_moranLocalIPerm <- function(gene,
                                 sfe,
                                 listw,
                                 nsim,
                                 alternative,
                                 na.action,
                                 zero.policy,
                                 mlvar,
                                 spChk,
                                 adjust.x,
                                 sample_Ei,
                                 iseed,
                                 no_repeat_in_row) {
  ## Select gene expression input
  x <- SummarizedExperiment::assay(sfe, "logcounts")[gene,]

  ## Check input validity
  .int_checkSAInput(x = x, listw = listw)

  ## Call the localmoran_perm function from spdep
  res <- spdep::localmoran_perm(x = x,
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

  ## Get attributes into a dataframe
  out <- data.frame(Ii = res[,"Ii"],
                    IiFDR = res[,dimnames(res)[[2]][5]],
                    IiClust = attr(res, "quadr")[,"pysal"]) # maybe give option to select between mean, median, and pysal

  return(out)
}


