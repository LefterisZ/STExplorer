gwrSTE <- function(gwr_method) {

}


#' Bandwidth selection for different GWR methods
#'
#' @description A function for automatic bandwidth selection to calibrate a
#' basic GWR model. It is a wrapper of different functions from the `GWmodel`
#' package, developed and maintained by Binbin Lu.
#'
#' @param formula Regression model formula of a \code{\link{formula}} object.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#' MetaSpatialFeatureExperiment.
#'
#' @param sample_id character string, TRUE or NULL specifying sample/image
#' identifier(s); here, TRUE is equivalent to all samples/images.
#'
#' @param approach Approach specified by \code{CV} for cross-validation
#' approach or by \code{AICc} corrected approach.
#'
#' @param kernel Function chosen as follows:
#'                \itemize{
#'                  \item{gaussian}{wgt = exp(-.5*(vdist/bw)^2);}
#'                  \item{exponential}{wgt = exp(-vdist/bw);}
#'                  \item{bisquare}{wgt = (1-(vdist/bw)^2)^2 if vdist < bw,
#' wgt=0 otherwise;}
#'                  \item{tricube}{wgt = (1-(vdist/bw)^3)^3 if vdist < bw,
#' wgt=0 otherwise;}
#'                  \item{boxcar}{wgt=1 if dist < bw, wgt=0 otherwise.}
#'                }
#'
#' @param adaptive If \code{TRUE}, calculate an adaptive kernel where the
#' bandwidth (\code{bw}) corresponds to the number of nearest neighbours
#' (i.e., adaptive distance); default is \code{FALSE}, where a fixed kernel is
#' found (bandwidth is a fixed distance).
#'
#' @param p The power of the Minkowski distance, default is 2, i.e., the
#' Euclidean distance.
#'
#' @param dMat A pre-specified distance matrix, it can be calculated by the
#' function \code{\link{gw.dist}}.
#'
#' @param family a description of the error distribution and link function to
#' be used in the model, which can be specified by “poisson” or “binomial”
#' description. Used only for the `gwr_method`: `"gwr-generalised"`
#' \code{\link[GWmodel]{bw.ggwr}}
#'
#' @param obs.tv a vector of time tags for each observation, which could be
#' numeric or of POSIXlt class. Used only for the `gwr_method`: `"gtwr"`
#' \code{\link[GWmodel]{bw.gtwr}}
#'
#' @param lamda an parameter between 0 and 1 for calculating spatio-temporal
#' distance. Used only for the `gwr_method`: `"gtwr"`
#' \code{\link[GWmodel]{bw.gtwr}}
#' @param lambda_lcr option for a globally-defined (constant) ridge parameter.
#' Default is lambda=0, which gives a basic GWR fit. Used only for
#' the `gwr_method`: `"gwr-lcr"` \code{\link[GWmodel]{bw.gtwr}}
#'
#' @param lambda.adjust a locally-varying ridge parameter. Default FALSE,
#' refers to: (i) a basic GWR without a local ridge adjustment (i.e. lambda=0,
#' everywhere); or (ii) a penalised GWR with a global ridge adjustment (i.e.
#' lambda is user-specified as some constant, other than 0 everywhere); if
#' TRUE, use cn.tresh to set the maximum condition number. For locations with a
#' condition number (for its local design matrix), above this user-specified
#' threshold, a local ridge parameter is found
#'
#' @param cn.thresh maximum value for condition number, commonly set between
#' 20 and 30
#'
#' @param t.units character string to define time unit. Used only for the
#' `gwr_method`: `"gtwr"`
#'
#' @param ksi a parameter between 0 and PI for calculating spatio-temporal
#' distance, see details in Wu et al. (2014). Used only for the
#' `gwr_method`: `"gtwr"`
#'
#' @param st.dMat a pre-specified spatio-temporal distance matrix. Used only
#' for the `gwr_method`: `"gtwr"`
#'
#' @param verbose logical variable to define whether show the selection
#' procedure. Used only for the `gwr_method`: `"gtwr"`
#'
#' @param parallel.method FALSE as default, and the calibration will be
#' conducted traditionally via the serial technique, "omp": multi-thread
#' technique with the OpenMP API, "cluster": multi-process technique with the
#' \pkg{parallel} package, "cuda": parallel computing technique with CUDA.
#'
#' @param parallel.arg If \code{parallel.method} is not FALSE, then set the
#' argument as follows: if \code{parallel.method} is "omp",
#' \code{parallel.arg} refers to the number of threads used, and its default
#' value is the number of cores - 1; if \code{parallel.method} is "cluster",
#' \code{parallel.arg} refers to the number of R sessions used, and its default
#' value is the number of cores - 1; if \code{parallel.method} is "cuda",
#' \code{parallel.arg} refers to the number of calibrations included in each
#' group, but note a too large value may cause the overflow of GPU memory.
#'
#' @note For a discontinuous kernel function, a bandwidth can be specified
#' either as a fixed (constant) distance or as a fixed (constant) number of
#' local data (i.e., an adaptive distance).  For a continuous kernel
#' function, a bandwidth can be specified either as a fixed distance or as a
#' 'fixed quantity that reflects local sample size'  (i.e., still an
#' 'adaptive' distance but the actual local sample size will be the sample
#' size as functions are continuous).  In practice, a fixed bandwidth suits
#' fairly regular sample configurations whilst an adaptive bandwidth suits
#' highly irregular sample configurations. Adaptive bandwidths ensure
#' sufficient (and constant) local information for each local calibration.
#' This note is applicable to all GW models.
#'
#' @return Returns the adaptive or fixed distance bandwidth.
#'
#' @importFrom GWmodel bw.ggwr bw.gtwr bw.gwr.lcr bw.gwr1
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @export

gwr_bwSTE <- function(gwr_method = c("basic", "gtwr",
                                     "gwr-lcr", "gwr-generalised"),
                      formula,
                      m_sfe,
                      sample_id,
                      approach = c("CV", "AIC"),
                      kernel = c("bisquare", "gaussian", "exponential",
                                 "tricube", "boxcar"),
                      family = c("poisson", "binomial"),
                      adaptive = FALSE,
                      p = 2,
                      dMat = NULL,
                      parallel.method = FALSE,
                      parallel.arg = NULL,
                      lamda = 0.05,
                      t.units = "auto",
                      ksi = 0,
                      st.dMat = NULL,
                      obs.tv = NULL,
                      verbose = TRUE,
                      lambda_lcr = 0,
                      lambda.adjust = FALSE,
                      cn.thresh = NA) {

  ## Check arguments
  method <- match.arg(gwr_method)

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Select dMat if not provided

  ## Prepare data for GWR

  ## Run GWR
  if (method == "basic") {
    return(GWmodel::bw.gwr(formula = formula,
                           data = data,
                           approach = approach,
                           kernel = kernel,
                           adaptive = adaptive,
                           p = p,
                           theta = 0,
                           longlat = FALSE,
                           dMat = dMat,
                           parallel.method = parallel.method,
                           parallel.arg = parallel.arg))

  } else if (method == "gtwr") {
    ## Automatic bandwidth selection to calibrate a GTWR model
    ## GTWR: Geographically and Temporally Weighted Regression
    return(GWmodel::bw.gtwr(formula = formula,
                            data = data,
                            obs.tv = obs.tv,
                            approach = approach,
                            kernel = kernel,
                            adaptive = adaptive,
                            p = p,
                            theta = 0,
                            longlat = FALSE,
                            lamda = lamda,
                            t.units = t.units,
                            ksi = ksi,
                            st.dMat = st.dMat,
                            verbose = verbose))

  } else if (method == "gwr-lcr") {
    return(GWmodel::bw.gwr.lcr(formula = formula,
                               data = data,
                               kernel = kernel,
                               lambda = lambda_lcr,
                               lambda.adjust = lambda.adjust,
                               cn.thresh = cn.thresh,
                               adaptive = adaptive,
                               p = p,
                               theta = 0,
                               longlat = FALSE,
                               dMat = dMat))

  } else if (method == "gwr-generalised") {
    ## Automatic bandwidth selection to calibrate a generalised GWR model
    return(GWmodel::bw.ggwr(formula = formula,
                            data = data,
                            family = family,
                            approach = approach,
                            kernel = kernel,
                            adaptive = adaptive,
                            family = family,
                            p = p,
                            theta = 0,
                            longlat = FALSE,
                            dMat = dMat))

  } else {
    stop("Invalid value for 'gwr_method'. Choose one of 'basic', 'gtwr',",
         " 'gwr-lcr', or 'gwr-generalised'.")
  }
}


