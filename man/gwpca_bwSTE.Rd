% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwpca_STE.R
\name{gwpca_bwSTE}
\alias{gwpca_bwSTE}
\title{Automatically find optimal bandwidth}
\usage{
gwpca_bwSTE(
  sfe,
  vars,
  k = 2,
  robust = FALSE,
  kernel = "gaussian",
  adaptive = FALSE,
  p = 2,
  dMat
)
}
\arguments{
\item{sfe}{A \code{SpatialFeatureExperiment} object.}

\item{vars}{A vector of variable names.}

\item{k}{The number of retained components. Default is 2.}

\item{robust}{If TRUE, apply robust GWPCA; otherwise, apply basic GWPCA.
Default is FALSE.}

\item{kernel}{The kernel function to use. Options include:
\itemize{
\item "gaussian": Gaussian kernel (default)
\item "exponential": Exponential kernel
\item "bisquare": Bisquare kernel
\item "tricube": Tricube kernel
\item "boxcar": Boxcar kernel
}}

\item{adaptive}{If TRUE, calculate an adaptive kernel where the bandwidth
corresponds to the number of nearest neighbors (adaptive distance). Default
is FALSE, where a fixed kernel is used (fixed distance).}

\item{p}{The order of the norm for Minkowski distance. If p = 1, it
represents Manhattan distance; if p = 2, it represents Euclidean distance.
Default is 2.}

\item{dMat}{The distance matrix. Default is NULL.}
}
\value{
Numeric, the optimal bandwidth.
}
\description{
A function to automatically find an optimal fixed or adaptive bandwidth.
This function is used to calibrate a basic or robust GWPCA via a
cross-validation approach.
}
\examples{
# Load a SpatialFeatureExperiment object
data(sfe)

# Prepare some parameters
vars <- rownames(sfe)[1:50]

# Find the optimal bandwidth using basic GWPCA
\dontrun{
bw <- gwpca_bwSTE(sfe = sfe, vars = vars)
bw
}

# Find the optimal bandwidth using robust GWPCA and an adaptive kernel
\dontrun{
bw <- gwpca_bwSTE(sfe = sfe, vars = vars,
robust = TRUE, adaptive = TRUE)
bw
}

# Find the optimal bandwidth using a pre-calculated distance matrix
\dontrun{
dMat <- gw.dist(dp.locat = spatialCoords(sfe), p = 2)
bw <- gwpca_bwSTE(sfe = sfe, vars = vars, dMat = dMat)
bw
}

}
\author{
Eleftherios (Lefteris) Zormpas
}
