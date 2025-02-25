% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_moranI.R
\name{moranGlobalI}
\alias{moranGlobalI}
\title{Compute Moran's I}
\usage{
moranGlobalI(
  m_sfe,
  sample_id = NULL,
  genes = TRUE,
  assay = "logcounts",
  n = NULL,
  S0 = NULL,
  zero.policy = NULL,
  NAOK = FALSE,
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

\item{n}{number of zones}

\item{S0}{global sum of weights}

\item{zero.policy}{Default is NULL. If not changed then internally, it is
set to \code{attr(listw, "zero.policy")} as set when \code{listw} was
created. If TRUE, assign zero to the lagged value of zones without
neighbours. If FALSE, assign NA. If attribute not set, use the global option
value.}

\item{NAOK}{if 'TRUE' pass 'NA' or 'NaN' or 'Inf' to foreign function}

\item{mc.cores}{Argument from \code{link[parallel]{mclapply}}. The number of
cores to use, i.e., at most how many child processes will be run
simultaneously. The option is initialized from environment variable MC_CORES
if set. Must be at least one, and parallelisation requires at least two
cores.}
}
\value{
An SFE object with the below two columns added in the \code{rowData}:
\item{moranI}{Moran's I}
\item{moranK}{sample kurtosis of x}
}
\description{
This function is a wrapper of the \code{\link[spdep]{moran}} function from
the \code{spdep} package developed by Roger Bivand. It is a function to compute
Moran's I, called by \code{moran.test} and \code{moran.mc}.
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

# Calculate Moran's I for this gene
moran <- moranGlobalI(x = gene_exp, listw = listw,zero.policy = TRUE)
}

}
\references{
Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 17.
}
\seealso{
\code{\link{moran.test}}, \code{\link{moran.mc}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
