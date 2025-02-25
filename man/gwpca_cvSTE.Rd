% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwpca_STE.R
\name{gwpca_cvSTE}
\alias{gwpca_cvSTE}
\title{Cross-validation contribution for GWPCA score}
\usage{
gwpca_cvSTE(
  bw,
  x,
  loc,
  dMat,
  k,
  robust = FALSE,
  kernel = "gaussian",
  adaptive = FALSE,
  cvContrib = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{bw}{Bandwidth for GWPCA.}

\item{x}{Data matrix.}

\item{loc}{Location matrix.}

\item{dMat}{Distance matrix.}

\item{k}{Number of neighbors for GWPCA.}

\item{robust}{Logical, indicating whether to use robust GWPCA.}

\item{kernel}{Kernel function for GWPCA.}

\item{adaptive}{Logical, indicating whether to use adaptive bandwidth.}

\item{cvContrib}{Logical, indicating whether to return a vector of
contributions for each location.}

\item{verbose}{Logical, indicating whether to display progress information.}
}
\value{
If cvContrib is FALSE, returns the sum of contributions to the score
statistic used in cross-validation. If cvContrib is TRUE, returns a vector of
contributions for each location. The cvContrib == FALSE is used when a
bandwidth is automatically selected using \code{gwpca_bwSTE} function. The
cvContrib == TRUE is used when we need a score for each location to identify
discrepencies and outliers.
}
\description{
Calculates the contribution of each observation to the score statistic used
in cross-validation for GWPCA.
}
\details{
The contribution of each observation to the score statistic used in
cross-validation for GWPCA is calculated. Outliers are taken to correspond to
high score (residual) values.
}
\examples{
# Example usage
# Set parameters
bw <- 0.5
k <- 5

# Create toy matrices
x <- matrix(rnorm(1000), ncol = 10)
loc <- abs(matrix(rnorm(200), ncol = 2))
dMat <- GWmodel::gw.dist(loc)

# Calculate
gwpca_cvSTE(bw, x, loc, dMat, k)

}
\seealso{
\code{\link{wpca}}, \code{\link{rwpca}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
