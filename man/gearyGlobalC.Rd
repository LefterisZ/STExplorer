% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_gearyC.R
\name{gearyGlobalC}
\alias{gearyGlobalC}
\title{Calculate Geary's C Global Statistic}
\usage{
gearyGlobalC(
  m_sfe,
  sample_id = NULL,
  genes = TRUE,
  assay = "logcounts",
  n = NULL,
  n1 = NULL,
  S0 = NULL,
  zero.policy = NULL,
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

\item{n}{Number of zones. If not provided, it is calculated as
length(listw$neighbours).}

\item{n1}{\code{n - 1}. If not provided, it is calculated as
\code{length(listw$neighbours) - 1}.}

\item{S0}{Global sum of weights. If not provided, it is calculated as
\code{spdep::Szero(listw)}.}

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
}
\value{
An SFE object with the below two columns added in the \code{rowData}:
\itemize{
\item gearyC_stat - the value of the observed Geary's C.
\item gearyK - sample kurtosis of gene.
}
}
\description{
This function is a wrapper of the \code{\link[spdep]{geary}} function from
the \code{spdep} package developed by Roger Bivand. It computes Geary's C
statistic for spatial autocorrelation using the spdep package. In its core
the function asks for two arguments: 1) \code{x} - A numeric vector of
observations. Must be the same length as the neighbours list in listw, and
2) \code{listw} A listw object created, for example, by
\code{link[spdep]{nb2listw}} or \code{link[spdep]{nb2listwdist}}.
}
\details{
If n, n1, and S0 are not provided, they are calculated based on the
\code{listw} object.
}
\examples{
\dontrun{
# Export an SFE object from an MSFE object
sfe <- getSFE(msfe, "JBO019")
# A vector of gene expression of one gene over all locations
gene_exp <- assay(sfe, "logcounts")[1,]
# A listw object containing distance-based spatial weights for neighbors
listw <- listw
# Calculate Geary's C for this gene
geary <- gearyGlobalC(x = gene_exp, listw = listw, zero.policy = TRUE)
}

}
\seealso{
\code{\link[spdep]{geary}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{autocorrelation}
\keyword{geary}
\keyword{listw}
\keyword{spatial}
