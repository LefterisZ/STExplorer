% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_moranI.R
\name{moranLocalI}
\alias{moranLocalI}
\alias{moranLocalIPerm}
\title{Local Moran's I statistic}
\usage{
moranLocalI(
  m_sfe,
  sample_id = NULL,
  genes = TRUE,
  assay = "logcounts",
  zero.policy = NULL,
  na.action = na.fail,
  conditional = TRUE,
  alternative = "two.sided",
  mlvar = TRUE,
  spChk = NULL,
  adjust.x = FALSE,
  mc.cores = getOption("mc.cores", 2L)
)

moranLocalIPerm(
  m_sfe,
  sample_id = NULL,
  genes = TRUE,
  assay = "logcounts",
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
  mc.cores = getOption("mc.cores", 2L)
)
}
\arguments{
\item{m_sfe}{An object of class SpatialFeatureExperiment or
MetaSpatialFeatureExperiment.}

\item{sample_id}{A character vector specifying the sample IDs to include
(only relevant if a MetaSpatialFeatureExperiment has been provided in the
\code{m_sfe} argument).}

\item{genes}{TRUE or a named character vector with gene names for which the
SA statistic needs to be calculated. If left to TRUE then the SA statistic
is calculated for every gene.}

\item{assay}{the counts assay to use. Defaults to "logcounts".}

\item{zero.policy}{Default is NULL. If not changed then internally, it is
set to \code{attr(listw, "zero.policy")} as set when \code{listw} was
created. If TRUE, assign zero to the lagged value of zones without
neighbours. If FALSE, assign NA. If attribute not set, use the global option
value.}

\item{na.action}{A function (default \code{na.fail}), can
also be \code{na.omit} or \code{na.exclude} - in these
cases the weights list will be subsetted to remove NAs
in the data. It may be necessary to set zero.policy to
TRUE because this subsetting may create no-neighbor
observations. Note that only weights lists created
without using the glist argument to \code{nb2listw} may
be subsetted. If \code{na.pass} is used, zero is
substituted for NA values in calculating the spatial
lag. (Note that na.exclude will only work properly
starting from R 1.9.0, na.omit and na.exclude assign the
wrong classes in 1.8.*)}

\item{conditional}{Default TRUE: expectation and variance
are calculated using the conditional randomization null
(Sokal 1998 Eqs. A7 & A8). Elaboration of these changes
available in Sauer et al. (2021). If FALSE: expectation
and variance are calculated using the total
randomization null (Sokal 1998 Eqs. A3 & A4).}

\item{alternative}{A character string specifying the
alternative hypothesis, must be one of greater, less or
two.sided (default).}

\item{mlvar}{Default TRUE: values of local Moran's I are
reported using the variance of the variable of interest
(sum of squared deviances over n), but can be reported
as the sample variance, dividing by (n-1) instead; both
are used in other implementations.}

\item{spChk}{Should the data vector names be checked
against the spatial objects for identity integrity,
TRUE, or FALSE, default NULL to use
\code{get.spChkOption()}.}

\item{adjust.x}{Default FALSE, if TRUE, x values of
observations with no neighbors are omitted in the mean
of x.}

\item{mc.cores}{Argument from \code{link[parallel]{mclapply}}. The number of
cores to use, i.e., at most how many child processes will be run
simultaneously. The option is initialized from environment variable MC_CORES
if set. Must be at least one, and parallelisation requires at least two
cores.}

\item{nsim}{Default 499, number of conditional
permutation simulations.}

\item{sample_Ei}{Default TRUE; if conditional permutation,
use the sample $E_i$ values, or the analytical values,
leaving only variances calculated by simulation.}

\item{iseed}{Default NULL, used to set the seed for
possible parallel RNGs.}

\item{no_repeat_in_row}{Default \code{FALSE}, if
\code{TRUE}, sample conditionally in each row without
replacements to avoid duplicate values,
\url{https://github.com/r-spatial/spdep/issues/124}.}
}
\value{
An SFE object with the results added in the \code{localResults} slot
which contains a DataFrame named \code{localMoranI}/\code{localMoranIPerm} that
contains the per location Ii statistics values for each gene:
\itemize{
\item{ENSG***.Ii}{Numeric vector containing Local Moran statistic
with attribute}
\item{ENSG***.IiFDR}{(Only for gearyLocalCPerm) \code{rank()} and
\code{punif()} of observed statistic rank for \link{0, 1} p-values using
\code{alternative=}}
\item{ENSG***.IiCluster}{Low/High cluster for each location. It can take
values: "Low-Low", "High-High", "Low-High", "High-Low"}
}
}
\description{
This function is a wrapper of the \code{\link[spdep]{moran}} function from
the \code{spdep} package developed by Roger Bivand. The local spatial statistic
Moran's I is calculated for each zone based on the spatial weights object
used. The values returned include a Z-value, and may be used as a diagnostic
tool.
}
\details{
The values of local Moran's I are divided by the variance (or
sample variance) of the variable of interest to accord with Table 1, p. 103,
and formula (12), p. 99, in Anselin (1995), rather than his formula (7), p.
98. The variance of the local Moran statistic is taken from Sokal et al.
(1998) p. 334, equations 4 & 5 or equations 7 & 8 located depending on user
specification. By default, the implementation divides by n, not (n-1) in
calculating the variance and higher moments. Conditional code contributed by
Jeff Sauer and Levi Wolf.
}
\note{
Conditional permutations added for comparative purposes; permutations
are over the whole data vector omitting the observation itself. For p-value
adjustment, use \code{p.adjust()} or \code{p.adjustSP()} on the output
vector.
}
\examples{
\dontrun{
# Example usage:

# Export an SFE object from an MSFE object
sfe <- getSFE(msfe, "JBO019")

# A vector of gene expression of one gene over all locations
gene_exp <- assay(sfe, "logcounts")[1,]

# A listw object containing distance-based spatial weights for neighbours
listw <- listw

# Calculate local Moran's I for this gene
moran <- moranLocalI(x = gene_exp, listw = listw,zero.policy = TRUE)
}

\dontrun{
# Calculate local Moran's I for this gene
moran <- moranLocalIPerm(x = gene_exp, listw = listw, nsim = 499,
zero.policy = TRUE)
}

}
\references{
Anselin, L. 1995. Local indicators of spatial association,
Geographical Analysis, 27, 93--115;
Getis, A. and Ord, J. K. 1996 Local spatial statistics: an overview. In P.
Longley and M. Batty (eds) \emph{Spatial analysis: modelling in a GIS
environment} (Cambridge: Geoinformation International), 261--277;
Sokal, R. R, Oden, N. L. and Thomson, B. A. 1998. Local Spatial
Autocorrelation in a Biological Model. Geographical Analysis, 30. 331--354;
Bivand RS, Wong DWS 2018 Comparing implementations of global and local
indicators of spatial association. TEST, 27(3), 716--748 \doi{10.1007/s11749-
018-0599-x};
Sauer, J., Oshan, T. M., Rey, S., & Wolf, L. J. 2021. The Importance of Null
Hypotheses: Understanding Differences in Local Moran’s under
Heteroskedasticity. Geographical Analysis. \doi{doi:10.1111/gean.12304}
Bivand, R. (2022), R Packages for Analyzing Spatial Data: A Comparative Case
Study with Areal Data. Geographical Analysis, 54(3), 488-518. \doi{10.1111/
gean.12319}
}
\seealso{
\code{\link{localG}}

\code{\link{localmoran_perm}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
