#' Geographically Weighted Regression (GWR)
#'
#' A function to perform Geographically Weighted Regression (GWR) using different
#' methods specified by the user.
#'
#' @param gwr_method Character vector specifying the GWR method to use. Possible
#'   values are: "basic" for basic GWR, "gtwr" for Geographically and Temporally
#'   Weighted Regression, "gwr-lcr" for Locally Compensated Ridge GWR (GWR-LCR),
#'   and "gwr-generalised" for Generalised GWR modeλ.
#'
#' @param formula Regression model formula of a \code{\link{formula}} object.
#'
#' @param m_sfe An object of class SpatialFeatureExperiment or
#'   MetaSpatialFeatureExperiment.
#'
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#'   identifier(s); here, TRUE is equivalent to all samples/images.
#'
#' @param bw Bandwidth for the GWR.
#'
#' @param kernel Function chosen as follows:
#'   \itemize{
#'     \item{gaussian}{wgt = exp(-.5*(vdist/bw)^2);}
#'     \item{exponential}{wgt = exp(-vdist/bw);}
#'     \item{bisquare}{wgt = (1-(vdist/bw)^2)^2 if vdist < bw,
#'       wgt=0 otherwise;}
#'     \item{tricube}{wgt = (1-(vdist/bw)^3)^3 if vdist < bw,
#'       wgt=0 otherwise;}
#'     \item{boxcar}{wgt=1 if dist < bw, wgt=0 otherwise.}
#'   }
#'
#' @param assay The counts assay to use. Defaults to "logcounts".
#'
#' @param adaptive If \code{TRUE}, calculate an adaptive kernel where the
#'   bandwidth (\code{bw}) corresponds to the number of nearest neighbours
#'   (i.e., adaptive distance); default is \code{FALSE}, where a fixed kernel
#'   is found (bandwidth is a fixed distance).
#'
#' @param p The power of the Minkowski distance, default is 2, i.e., the
#'   Euclidean distance.
#'
#' @param dMat A pre-specified distance matrix, it can be calculated and added
#'   to the SFE object by the \code{\link{addDistMat}} function. If `NULL`, it
#'   fetches the first (or only) distance matrix from the SFE's metadata slot.
#'   If you have multiple distance matrices then you can use one of the
#'   `"euclidean"`, `"manhattan"`, or `"minkowski"` to select the one you
#'   prefer. Defaults to `NULL`. If you want to calculate a different distance
#'   matrix, you can do so by using the \code{\link[GWmodel]{gw.dist}} function.
#'   Then, provide the matrix as a value to this argument.
#'
#' @param F123.test If TRUE, conduct three separate F-tests according to
#'   Leung et al. (2000).
#'
#' @param cv If TRUE, cross-validation data will be calculated and returned in
#'   the output Spatial*DataFrame.
#'
#' @param W.vect Default NULL, if given it will be used to weight the distance
#'   weighting matrix.
#'
#' @param parallel.method FALSE as default, and the calibration will be
#'   conducted traditionally via the serial technique, "omp": multi-thread
#'   technique with the OpenMP API, "cluster": multi-process technique with
#'   the \pkg{parallel} package, "cuda": parallel computing technique with CUDA.
#'
#' @param parallel.arg If \code{parallel.method} is not FALSE, then set the
#'   argument as follows: if \code{parallel.method} is "omp",
#'   \code{parallel.arg} refers to the number of threads used, and its default
#'   value is the number of cores - 1; if \code{parallel.method} is "cluster",
#'   \code{parallel.arg} refers to the number of R sessions used, and its
#'   default value is the number of cores - 1; if \code{parallel.method} is
#'   "cuda", \code{parallel.arg} refers to the number of calibrations included
#'   in each group, but note a too large value may cause the overflow of GPU
#'   memory.
#'
#' @param lamda A parameter between 0 and 1 for calculating spatio-temporal
#'   distance. Used only for the `gwr_method`: `"gtwr"`.
#'
#' @param t.units Character string to define time unit. Used only for the
#'   `gwr_method`: `"gtwr"`.
#'
#' @param ksi A parameter between 0 and PI for calculating spatio-temporal
#'   distance, see details in Wu et al. (2014). Used only for the
#'   `gwr_method`: `"gtwr"`.
#'
#' @param st.dMat A pre-specified spatio-temporal distance matrix. Used only
#'   for the `gwr_method`: `"gtwr"`.
#'
#' @param obs.tv A vector of time tags for each observation, which could be
#'   numeric or of POSIXlt class. Used only for the `gwr_method`: `"gtwr"`.
#'
#' @param lambda_lcr Option for a globally-defined (constant) ridge parameter.
#'   Default is lambda=0, which gives a basic GWR fit. Used only for
#'   the `gwr_method`: `"gwr-lcr"`.
#'
#' @param lambda.adjust A locally-varying ridge parameter. Default FALSE,
#'   refers to: (i) a basic GWR without a local ridge adjustment (i.e. lambda=0,
#'   everywhere); or (ii) a penalised GWR with a global ridge adjustment (i.e.
#'   lambda is user-specified as some constant, other than 0 everywhere); if
#'   TRUE, use cn.tresh to set the maximum condition number. For locations with
#'   a condition number (for its local design matrix), above this
#'   user-specified threshold, a local ridge parameter is found.
#'
#' @param cn.thresh Maximum value for condition number, commonly set between
#'   20 and 30.
#'
#' @return Returns the GWR model object.
#' @author Eleftherios Zormpas
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' gwrSTE(gwr_method = "basic",
#'        formula = Y ~ X1 + X2,
#'        m_sfe = sfe_object,
#'        sample_id = TRUE,
#'        bw = 100,
#'        kernel = "gaussian",
#'        adaptive = FALSE,
#'        p = 2,
#'        F123.test = FALSE,
#'        cv = TRUE)
#' }
#' @export
gwrSTE <- function(gwr_method = c("basic", "gtwr",
                                  "gwr-lcr", "gwr-generalised"),
                   formula,
                   m_sfe,
                   sample_id,
                   bw,
                   kernel = c("bisquare", "gaussian", "exponential",
                              "tricube", "boxcar"),
                   assay = "logcounts",
                   adaptive = FALSE,
                   p = 2,
                   dMat = NULL,
                   F123.test = FALSE,
                   cv = TRUE,
                   W.vect = NULL,
                   parallel.method = FALSE,
                   parallel.arg = NULL,
                   lamda = 0.05,
                   t.units = "auto",
                   ksi = 0,
                   st.dMat = NULL,
                   obs.tv = NULL,
                   lambda_lcr = 0,
                   lambda.adjust = FALSE,
                   cn.thresh = NA,
                   scale = TRUE) {
  ## Check arguments
  method <- match.arg(gwr_method)

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Select dMat if not provided
  ## Call internal function to select and obtain the distance matrix
  dMat <- .int_checkDMat(dMat, sfe)

  ## Prepare data for GWR
  data <- .int_getGWRdata(sfe, formula = formula, type = "hex", assay = assay)

  ## Scale data for Regression
  if (scale) {
    for (i in 1:ncol(data)) {
      data[[i]] <- base::scale(data[[i]], center = TRUE, scale = TRUE)
    }
  }

  ## Run GWR
  if (method == "basic") {
    ## basic GWR
    res <- GWmodel::gwr.basic(formula = formula,
                              data = data,
                              # regression.points,
                              bw = bw,
                              kernel = kernel,
                              adaptive = adaptive,
                              p = p,
                              theta = 0,
                              longlat = FALSE,
                              dMat = dMat,
                              F123.test = F123.test,
                              cv = cv,
                              W.vect = W.vect,
                              parallel.method = parallel.method,
                              parallel.arg = parallel.arg)

  } else if (method == "gtwr") {
    ## GTWR: Geographically and Temporally Weighted Regression
    res <- GWmodel::gtwr(formula = formula,
                         data = data,
                         # regression.points,
                         obs.tv = obs.tv,
                         # reg.tv,
                         st.bw = bw,
                         kernel = kernel,
                         adaptive = adaptive,
                         p = p,
                         theta = 0,
                         longlat = FALSE,
                         lamda = lamda,
                         t.units = t.units,
                         ksi = ksi,
                         st.dMat = st.dMat)

  } else if (method == "gwr-lcr") {
    ## Locally compensated ridge GWR (GWR-LCR)
    res <- GWmodel::gwr.lcr(formula = formula,
                            data = data,
                            # regression.points,
                            bw = bw,
                            kernel = kernel,
                            lambda = lambda_lcr,
                            lambda.adjust = lambda.adjust,
                            cn.thresh = cn.thresh,
                            adaptive = adaptive,
                            p = p,
                            theta = 0,
                            longlat = FALSE,
                            cv = cv,
                            dMat = dMat)

  } else if (method == "gwr-generalised") {
    ## Generalised GWR model
    res <- GWmodel::ggwr.basic(formula = formula,
                               data = data,
                               family = family,
                               # regression.points,
                               bw = bw,
                               kernel = kernel,
                               adaptive = adaptive,
                               cv = cv,
                               tol = tol,
                               maxiter = maxiter,
                               p = p,
                               theta = 0,
                               longlat = FALSE,
                               dMat = dMat,
                               dMat1 = NULL)

  }

  ## The lm function used internally for the global regression stores also an
  ## attribute named `.Environment` which can be very big in size taking up a
  ## lot of memory. Also stores technically the exact same information in two
  ## slots, one named `terms` and one named `model`. Here we remove the `terms`
  ## slot completely and we remove the .Environment` attribute from the `model`
  ## slot.
  res$lm$terms <- NULL
  attr(attr(res[["lm"]][["model"]], "terms"), ".Environment") <- NULL
  res$lm$terms <- attr(res[["lm"]][["model"]], "terms")

  ## Return the object
  return(res)
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
#' @param dMat A pre-specified distance matrix, it can be calculated and added
#' to the SFE object by the \code{\link{addDistMat}} function. If `NULL`, it
#' fetches the first (or only) distance matrix from the SFE's metadata slot.
#' If you have multiple distance matrices then you can use one of the
#' `"euclidean"`, `"manhattan"`, or `"minkowski"` to select the one you prefer.
#' Defaults to `NULL`. If you want to calculate a different distance matrix,
#' you can do so by using the \code{\link[GWmodel]{gw.dist}} function. Then,
#' provide the matrix as a value to this argument.
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
                      adaptive = FALSE,
                      p = 2,
                      dMat = NULL,
                      parallel.method = FALSE,
                      parallel.arg = NULL,
                      family = c("poisson", "binomial"),
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
  ## Call internal function to select and obtain the distance matrix
  dMat <- .int_checkDMat(dMat, sfe)

  ## Prepare data for GWR
  data <- .int_getGWRdata(sfe, type = "hex")

  ## Run GWR badnwidth selection
  if (method == "basic") {
    ## Bandwidth selection for basic GWR
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
    ## Bandwidth selection for locally compensated ridge GWR (GWR-LCR)
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
                            p = p,
                            theta = 0,
                            longlat = FALSE,
                            dMat = dMat))

  }
}


#' Extract SDF from gwr object
#'
#' This function extracts the SDF data frame from within a gwr class object and
#' converts it to class sf.
#'
#' @param gwr A GWR object as created by `gwrSTE`.
#'
#' @return A data frame of class sf.
#'
#' @importFrom sf st_as_sf
#'
#' @rdname gwr_toSF
#' @author Eleftherios (Lefteris) Zormpas
#' @export
gwr_toSF <- function(gwr) {
  sf::st_as_sf(gwr$SDF)
}



#' GWR summary table
#'
#' This function generates a summary table for a Geographically Weighted
#' Regression (GWR) object.
#'
#' @param gwr A GWR object as created by `gwrSTE`, or the SDF data frame from
#'            within the GWR object converted to an SF class object using
#'            `gwr_toSF`.
#' @param stat Character string. It must be a column name from the `gwr$SDF`
#'             data frame. It is essentially, the statistic for which to
#'             provide a summary.
#'
#' @return A summary table containing descriptive statistics for the GWR object.
#'
#' @rdname gwr_stats
#' @author Eleftherios Zormpas
#' @export
gwr_stats <- function(gwr, stat = NULL) {
  UseMethod("gwr_stats")
}

#' @rdname gwr_stats
#' @export
gwr_stats.gwrm <- function(gwr) {
  n <- .int_countElements(deparse1(gwr$lm$terms[[3]]))
  gwr.tab <- apply(gwr$SDF@data[, 1:(5 + n)], 2, summary)
  gwr.tab <- round(gwr.tab, 1)
  gwr.tab <- t(gwr.tab[,1:(1 + n)])
  return(gwr.tab)
}

#' @rdname gwr_stats
#' @export
gwr_stats.SF <- function(gwr, stat) {
  ## Summary statistics for Local_R2
  summary_stats <- summary(gwr[[stat]])
  cat("Summary Statistics for Local R-squared Values:\n")
  return(summary_stats)
}


# ---------------------------------------------------------------------------- #
#  ################# INTERNAL FUNCTIONS ASSOCIATED WITH GWR #################
# ---------------------------------------------------------------------------- #
#' Internal: Get the provided distance matrix
#'
#' This function checks if a distance matrix (\code{dMat}) is provided. If
#' provided, it returns the distance matrix; otherwise, it returns \code{NULL}.
#'
#' @param dMat A matrix representing the distance matrix.
#' @return A distance matrix if provided, otherwise \code{NULL}.
#' @keywords internal
.int_getProvidedDMat <- function(dMat) {
  if (is.matrix(dMat)) {
    return(dMat)
  }
  NULL  # Return NULL if dMat is not provided
}

#' Internal: Get the distance matrix from spatial feature object's metadata
#'
#' This function attempts to obtain the distance matrix from the metadata of a
#' spatial feature experiment object (\code{sfe}) when the dMat argument
#' provided in the exported function is left to NULL. If so, it returns the
#' first distance matrix; otherwise, it returns \code{NULL}.
#'
#' @param sfe The spatial feature experiment object.
#' @param dMat expected to be NULL
#' @return A distance matrix if found in metadata, otherwise \code{NULL}.
#' @keywords internal
.int_getMetadataDMat <- function(dMat, sfe) {
  if (is.null(dMat)) {
    return(metadata(sfe)[["dMat"]][[1]])
  }
  NULL  # Return NULL if dMat is not present in metadata
}

#' Internal: Get the distance matrix based on the specified type
#'
#' This function attempts to obtain the distance matrix from the metadata of a
#' spatial feature object (\code{sfe}) based on the specified type. If found,
#' it returns the distance matrix; otherwise, it returns \code{NULL}.
#'
#' @param sfe The spatial feature experiment object.
#' @param type The type of distance matrix ('euclidean', 'manhattan', or
#' 'minkowski').
#' @return A distance matrix if found based on the specified type, otherwise
#' \code{NULL}.
#' @keywords internal
.int_getTypedDMat <- function(sfe, type) {
  if (type %in% c("euclidean", "manhattan")) {
    return(metadata(sfe)[["dMat"]][[type]])
  } else if (type == "minkowski") {
    mwski <- grep("minkowski", names(metadata(sfe)[["dMat"]]), value = TRUE)
    return(metadata(sfe)[["dMat"]][[mwski]])
  }
  NULL  # Return NULL if type is not recognized
}

#' Internal: select and obtain the distance matrix
#'
#' This function selects and obtains the distance matrix (\code{dMat}) to be
#' used. It first checks if dMat is provided, then tries to obtain it from the
#' metadata of a spatial feature object (\code{sfe}), and finally attempts to
#' obtain it based on the specified type. If all attempts fail, it raises an
#' error asking the user to add a distance matrix using the \code{addDistMat}
#' function.
#'
#' @param dMat A matrix representing the distance matrix, or a string
#' indicating the type of distance matrix ('euclidean', 'manhattan', or
#' 'minkowski').
#' @param sfe The spatial feature experiment object from which to obtain the
#' distance matrix if not provided explicitly.
#' @return A distance matrix.
#' @keywords internal
.int_checkDMat <- function(dMat, sfe) {
  ## Check if dMat is provided
  dMatrix <- .int_getProvidedDMat(dMat = dMat)
  if (!is.null(dMatrix)) {
    return(dMatrix)
  }

  ## Try obtaining from metadata
  dMatrix <- .int_getMetadataDMat(dMat = dMat, sfe = sfe)
  if (!is.null(dMatrix)) {
    return(dMatrix)
  }

  ## Try obtaining based on type
  dMatrix <- .int_getTypedDMat(sfe = sfe, type = dMat)
  if (!is.null(dMatrix)) {
    return(dMatrix)
  }

  ## If all fails, raise an error
  stop("dMat argument is left to NULL and no dMat is present in the ",
       "SFE's metadata slot. Please use the `addDistMat` function to add ",
       "a distance matrix into the SFE objects.")
}


#' Internal: Prepare Data for Geographically Weighted Regression (GWR)
#'
#' An internal function to prepare data for Geographically Weighted Regression
#' (GWR) by retrieving geometries and converting data to the appropriate format.
#'
#' @param sfe An object of class SpatialFeatureExperiment.
#' @param type A character string specifying the type of geometry to use for GWR.
#' Options are "spot" (spot geometry), "hex" (hexagon geometry), or "centroid"
#' (centroid geometry).
#' @param formula The formula for the regression.
#' @param assay the counts assay to use
#'
#' @return A SpatialPointsDataFrame or SpatialPolygonsDataFrame object
#' containing the prepared data for GWR.
#'
#' @details This function prepares the data for Geographically Weighted Regression
#' (GWR) by retrieving the specified type of geometry from a SpatialFeatureExperiment
#' object and converting the data to the appropriate format. The available options
#' for the type parameter are "spot" (spot geometry), "hex" (hexagon geometry),
#' and "centroid" (centroid geometry).
#'
#' @author Eleftherios Zormpas
#' @rdname dot-int_getGWRdata
#' @importFrom sf st_as_sf
#' @importFrom dplyr left_join
#'
.int_getGWRdata <- function(sfe, formula, type, assay) {
  ## Get geometries
  polygons <- data.frame(geometry = .int_selectGeom(sfe = sfe, type = type),
                         Barcode = rownames(colGeometries(sfe)[[1]]))

  ## Extract values from formula
  ## Split by "+" or "~" and remove whitespace
  vars <- strsplit(formula, "\\s*(\\+\\s*|\\~\\s*)\\s*")[[1]]
  var_gs <- grepl("ENS", vars)
  var_loc <- vars %in% colnames(colData(sfe))
  var_non_gene_assay <- !var_gs & !var_loc

  ## Check whether the formula contains EnsgIDs, loc-specific values or both
  if (sum(var_gs) > 0 && sum(var_loc) == 0 && sum(var_non_gene_assay) == 0) {
    dt <- .int_getGWRdataGene(sfe = sfe, variables = vars[var_gs], assay = assay)
  } else if (sum(var_loc) > 0 && sum(var_gs) == 0 && sum(var_non_gene_assay) == 0) {
    dt <- .int_getGWRdataLoc(sfe = sfe, variables = vars[var_loc])
  } else if (sum(var_non_gene_assay) > 0 && sum(var_gs) == 0 && sum(var_loc) == 0) {
    dt <- .int_getGWRdataNonGeneAssay(sfe = sfe, variables = vars[var_non_gene_assay], assay = assay)
  } else {
    ## Get data for gene (if present)
    dt_gene <- if (sum(var_gs) > 0)
      .int_getGWRdataGene(sfe = sfe,
                          variables = vars[var_gs],
                          assay = assay) else NULL
    ## Get data for location (if present)
    dt_loc <- if (sum(var_loc) > 0)
      .int_getGWRdataLoc(sfe = sfe,
                         variables = vars[var_loc]) else NULL
    ## Get data for non-gene variables (if present)
    dt_non_gene_assay <- if (sum(var_non_gene_assay) > 0)
      .int_getGWRdataNonGeneAssay(sfe = sfe,
                                  variables = vars[var_non_gene_assay],
                                  assay = assay) else NULL

    ## Filter out NULL values before using Reduce
    dfs <- list(dt_gene, dt_loc, dt_non_gene_assay)
    non_null_dfs <- Filter(Negate(is.null), dfs)

    dt <- Reduce(function(x, y) dplyr::left_join(x, y, by = "Barcode"),
                 non_null_dfs)
  }

  ## Merge with polygons and convert to sp format
  dt <- dt %>%            # barcodes from row names to column
    dplyr::left_join(polygons, by = "Barcode") %>% # merge with geometries
    tibble::column_to_rownames(var = "Barcode") %>% # barcodes back to row names
    sf::st_as_sf(., sf_column_name = "geometry") %>% # data frame to sf
    as(., "Spatial")

  return(dt)
}


#' Internal: Retrieve Gene Data for GWR
#'
#' An internal function to extract gene expression data for use in
#' Geographically Weighted Regression (GWR).
#'
#' @param sfe An object of class SpatialFeatureExperiment.
#' @param variables A character vector of gene identifiers (e.g., ENSG IDs) to
#' be included in the GWR.
#' @param assay A character string specifying the assay to use for gene
#' expression data.
#'
#' @return A data frame with barcodes as row names and the specified gene
#' expression values as columns.
#'
#' @details This function retrieves gene expression data from the specified
#' assay of a SpatialFeatureExperiment object, selecting only the genes
#' specified in the variables parameter. The data is transposed to have
#' barcodes as row names and gene expressions as columns, suitable for further
#' processing in GWR.
#'
#' @author Eleftherios Zormpas
#' @rdname dot-int_getGWRdataGene
#' @importFrom dplyr select
#' @importFrom tibble rownames_to_column
#'
.int_getGWRdataGene <- function(sfe, variables, assay) {
  ## Fetch gene counts
  dt <- assay(sfe, assay) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::select(variables) %>% # select the genes in the formula
    tibble::rownames_to_column(var = "Barcode") # barcodes from row names to column

  return(dt)
}


#' Internal: Retrieve Location Data for GWR
#'
#' An internal function to extract location-specific variables for use in
#' Geographically Weighted Regression (GWR).
#'
#' @param sfe An object of class SpatialFeatureExperiment.
#' @param variables A character vector of column names from colData to be
#' included in the GWR.
#'
#' @return A data frame with barcodes as row names and the specified location
#' variables as columns.
#'
#' @details This function extracts location-specific data from the colData of
#' a SpatialFeatureExperiment object. It selects only the variables specified
#' in the variables parameter. The resulting data frame is formatted with
#' barcodes as row names, ready for integration into GWR analysis.
#'
#' @author Eleftherios Zormpas
#' @rdname dot-int_getGWRdataLoc
#' @importFrom dplyr select
#' @importFrom tibble rownames_to_column
#'
.int_getGWRdataLoc <- function(sfe, variables) {
  ## Fetch location variables
  dt <- colData(sfe) %>%
    as.data.frame() %>%
    dplyr::select(variables) %>% # select the location vars from the formula
    tibble::rownames_to_column(var = "Barcode") # sf to sp - GWmodel still works with sp objects

  return(dt)
}


#' Internal: Retrieve Non-Gene Assay Data for GWR
#'
#' An internal function to extract non-gene assay data for use in
#' Geographically Weighted Regression (GWR).
#'
#' @param sfe An object of class SpatialFeatureExperiment.
#' @param variables A character vector of non-gene assay identifiers to
#' be included in the GWR.
#' @param assay A character string specifying the assay to use for non-gene
#' data.
#'
#' @return A data frame with barcodes as row names and the specified non-gene
#' assay values as columns.
#'
#' @details This function retrieves non-gene data from the specified
#' assay of a SpatialFeatureExperiment object, selecting only the variables
#' specified in the variables parameter. The data is transposed to have
#' barcodes as row names and non-gene assay values as columns, suitable for further
#' processing in GWR.
#'
#' @importFrom dplyr select
#' @importFrom tibble rownames_to_column
#'
.int_getGWRdataNonGeneAssay <- function(sfe, variables, assay) {
  ## Fetch non-gene assay data
  dt <- assay(sfe, assay) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::select(variables) %>% # select the non-gene variables in the formula
    tibble::rownames_to_column(var = "Barcode") # barcodes from row names to column

  return(dt)
}


#' Internal: Function to count elements in each string
#'
#' This internal function splits a character string by the "+" symbol and counts
#' the number of elements in the resulting list.
#'
#' @param string A character string containing elements separated by "+".
#'
#' @return The number of elements in the string.
#'
#' @rdname dot-int_countElements
#' @author Eleftherios Zormpas
.int_countElements <- function(string) {
  ## Split by "+" and remove whitespace
  elements <- strsplit(string, "\\s*\\+\\s*")[[1]]
  return(length(elements))
}
