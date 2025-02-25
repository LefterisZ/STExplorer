% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_gearyC.R
\name{gearyLocalC}
\alias{gearyLocalC}
\alias{gearyLocalCPerm}
\title{Local Geary's C Statistic for Spatial Autocorrelation}
\usage{
gearyLocalC(
  m_sfe,
  sample_id = NULL,
  genes = TRUE,
  assay = "logcounts",
  ...,
  zero.policy = NULL,
  mc.cores = getOption("mc.cores", 2L)
)

gearyLocalCPerm(
  m_sfe,
  sample_id = NULL,
  genes = TRUE,
  assay = "logcounts",
  nsim = 999,
  alternative = "two.sided",
  ...,
  zero.policy = NULL,
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

\item{...}{other arguments passed to methods. Have a look at
\code{\link[spdep]{localC}}}

\item{zero.policy}{Default is NULL. If not changed then internally, it is
set to \code{attr(listw, "zero.policy")} as set when \code{listw} was
created. If TRUE, assign zero to the lagged value of zones without
neighbours. If FALSE, assign NA. If attribute not set, use the global option
value.}

\item{mc.cores}{Argument from \code{link[parallel]{mclapply}}. The number of
cores to use, i.e., at most how many child processes will be run
simultaneously. The option is initialized from environment variable MC_CORES
if set. Must be at least one, and parallelisation requires at least two
cores.}

\item{nsim}{The number of simulations to be used for permutation test.}

\item{alternative}{A character defining the alternative hypothesis. Must be
one of "two.sided", "less" or "greater".}

\item{iseed}{default NULL, used to set the seed for possible parallel RNGs.}

\item{no_repeat_in_row}{default FALSE, if TRUE, sample conditionally in each
row without replacements to avoid duplicate values,
https://github.com/r-spatial/spdep/issues/124}
}
\value{
An SFE object with the results added in the \code{localResults} slot
which contains a DataFrame named \code{localGearyC}/\code{localGearyCPerm} that
contains the per location Ci statistics values for each gene:
\itemize{
\item Ci or ENSG***.Ci - Numeric vector containing Local Geary statistic
with attribute
\item ENSG***.CiFDR - (Only for gearyLocalCPerm) \code{rank()} and
\code{punif()} of observed statistic rank for \link{0, 1} p-values using
\code{alternative=}
\item ENSG***.CiCluster - (Only for gearyLocalCPerm) Low/High cluster for
each location. It can take values: "Low-Low", "High-High", "Other Positive",
"Negative"
}

description
}
\description{
The Local Geary is a local adaptation of Geary's C
statistic of spatial autocorrelation. It uses squared
differences to measure dissimilarity, unlike the Local
Moran. Low values indicate positive spatial autocorrelation,
and large values refer to negative autocorrelation.
}
\details{
Inference for the Local Geary is based on permutation
approach, comparing the observed value to the reference
distribution under spatial randomness. \code{localC_perm()}
returns a pseudo p-value. This is not an analytical p-value
and should be used with care, based on permutations.

For more information about the usage cases and the
rest of the arguments visit the help page of \code{\link[spdep]{localC}}
Local Geary can be extended to a multivariate
context. When \code{x} is numeric, univariate Local Geary
is calculated. For multivariate, provide a list or matrix.

While not required in univariate context, standardized
Local Geary is calculated. Multivariate Local Geary is
\emph{always} standardized.

The univariate Local Geary is calculated as \eqn{c_i = \sum_j
w_{ij}(x_i - x_j)^2} and the multivariate Local Geary is
calculated as \eqn{c_{k,i} = \sum_{v=1}^{k} c_{v,i}} as
described in Anselin (2019).
}
\examples{
\dontrun{
# Export an SFE object from an MSFE object
sfe <- getSFE(msfe, "JBO019")
# A vector of gene expression of one gene over all locations
gene_exp <- assay(sfe, "logcounts")[1,]
# A listw object containing distance-based spatial weights for neighbours
listw <- listw
# Calculate Geary's C for this gene
geary <- gearyLocalC(x = gene_exp, listw = listw, zero.policy = TRUE)
}

\dontrun{
# Calculate  Geary's C and statistical significance for this gene
geary <- gearyLocalCPerm(x = gene_exp, listw = listw, nsim = 499,
zero.policy = TRUE)
}
}
\references{
Anselin, L. (2019). Local indicators of spatial
association—LISA. Spatial Statistics, 29, 10-27.
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{autocorrelation}
\keyword{geary}
\keyword{listw}
\keyword{spatial}
