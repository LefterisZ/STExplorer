#' Run GWPCA
#'
#' @name gwpcaSTE
#'
#' @description
#' A Geographically Weighted Principal Components Analysis function based on
#' \code{gwpca} function from \code{GWmodel} package. The function is re-written
#' and optimised to work faster and with a \code{SpatialFeatureExperimetn} (SFE)
#' object.
#'
#' @param sfe a \code{SpatialFeatureExperiment} (SFE) object.
#'
#' @param assay the counts assay to be used. Defaults to "logcounts".
#'
#' @param elocat a two-column numeric DataFrame object for providing evaluation
#' locations.
#'
#' @param vars a vector of variable names to be evaluated.
#'
#' @param p the order of the norm for Minkowski distance. If p = 1, represents
#' Manhattan distance; if p = 2, represents Euclidean distance.
#'
#' @param k the number of retained components; k must be less than the number of
#' variables
#'
#' @param bw bandwidth used in the weighting function, possibly calculated by
#' bw.gwpca;fixed (distance) or adaptive bandwidth(number of nearest neighbours)
#'
#' @param kernel function chosen as follows:
#' \itemize{
#'    \item gaussian: wgt = exp(-.5*(vdist/bw)^2)
#'    \item exponential: wgt = exp(-vdist/bw)
#'    \item bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise
#'    \item tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise
#'    \item boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' }
#'
#' @param adaptive if TRUE calculate an adaptive kernel where the bandwidth
#' corresponds to the number of nearest neighbours (i.e. adaptive distance);
#' default is FALSE, where a fixed kernel is found (bandwidth is a fixed
#' distance)
#'
#' @param scores if scores = TRUE, the scores of the supplied data on the
#' principal components will be calculated.
#'
#' @param robust if TRUE, robust GWPCA will be applied; otherwise basic GWPCA
#' will be applied.
#'
#' @param cv If TRUE, cross-validation data will be found that are used to
#' calculate the cross-validation score for the specified bandwidth.
#'
#' @param future Defaults to FALSE. If you already have set a \code{future}
#' backend then set it to true. It works as a switch to not confuse futures.
#'
#' @param strategy the future plan strategy to set. More info at future::plan.
#' It is used only if \code{future} argument == FALSE.
#'
#' @param workers the number of cores to be used. More info at future::plan.
#' It is used only if \code{future} argument == FALSE.
#'
#' @param verbose default TRUE. Show progress bar.
#'
#' @returns A list of class `gwpca`
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Load data -----------------------------
#' data(sfe)
#'
#' # Set parameters ------------------------
#' vars <- rownames(sfe)[1:50]
#' bw <- 650
#'
#' # Set special parameters ----------------
#' # Because it can be slow, we will select only a few locations to run over for
#' # this example.
#' elocat <-  colnames(sfe)[1:50]
#'
#' # Run GWPCA -----------------------------
#' pcagw <- gwpcaSTE(sfe = sfe, elocat = elocat, vars = vars, bw = bw)
#'
#' print(pcagw)
#'
#' @export
gwpcaSTE <- function(sfe,
                      assay = "logcounts",
                      elocat = NULL,
                      vars,
                      p = 2,
                      k = 20,
                      bw,
                      kernel = "gaussian",
                      adaptive = FALSE,
                      scores = FALSE,
                      robust = FALSE,
                      cv = TRUE,
                      future = FALSE,
                      strategy = "sequential",
                      workers = 1,
                      verbose = TRUE){

  ## Set handlers for a verbose approach
  if (verbose) {
    set.verbose(.global = TRUE,
                .handlers = c("progress"))
  }

  ## Set parallelisation if an appropriate strategy is supplied
  if (!future) {
    ## Fetch the user's current backend
    oplan <- plan()
    ## Set parallel or sequential backend
    set.parallel(.strategy = strategy, .workers = workers)
    ## doFuture is using a parallel-safe RNG method as per this:
    ## https://github.com/tidymodels/tune/issues/377 therefore we silence the
    ## warning.
    if (getDoParWorkers() > 1) {
      rlang::local_options(doFuture.rng.onMisuse = "ignore")
    }
  }
  ## Prepare the user for long waiting times
  locations <- nrow(spatialCoords(sfe))
  message("Locations to iterate over: ", locations)
  message("Please be patient. This will take a while...")
  mins <- round((locations * 2) / 60, 3)
  message("E.T.A. is approximatelly: ", mins, " minutes")
  message("NOTE: The time estimate is a function of the number of variables \n",
          "you provide (vars), the number of components you wish to keep (k),",
          "\nthe number of locations (`ncol(sfe)` or elocat) you have, \n",
          "whether you asked PCA scores to be calculated (scores = TRUE), \n",
          "the number ofcores you provided and the future plan you selected.\n",
          "---> Having these in mind, the E.T.A. estimate might be way off,\n",
          "     by overestimating (or underestimating) the time needed!")
  message("Sit back and relax!! :)")
  ## fetch start time
  s <- Sys.time()

  ## RUN GWPCA
  gwpca_result <- int.gwpca(.sfe = sfe,
                            .assay = assay,
                            .elocat = elocat,
                            .vars = vars,
                            .p = p,
                            .k = k,
                            .bw = bw,
                            .kernel = kernel,
                            .adaptive = adaptive,
                            .scores = scores,
                            .robust = robust,
                            .cv = cv,
                            .verbose = verbose)
  # print(gwpca_result)
  ## fetch end time
  e <- Sys.time()
  ## return difference
  message("Running GWPCA is done!")
  message("Time elapsed: ", round(difftime(e, s, units = "mins"), 3), " mins")


  # Reset the user's backend
  ## A future plan should never be called inside a function. But, because here
  ## we need to do so, we first save the user's backend and then we set it back.
  if (!future) {
    ## Stop cluster
    if (strategy == "cluster" |
        strategy == "clustAuto" |
        inherits(strategy, "ClusterFuture")) {
      parallel::stopCluster(workers)
    }
    ## Reset plan
    plan(oplan)
    plan()
  }

  ## Return the list
  return(gwpca_result)
}

# ---------------------------------------------------------------------------- #
#' Cross-validation contribution for GWPCA score
#'
#' Calculates the contribution of each observation to the score statistic used
#' in cross-validation for GWPCA.
#'
#' @name gwpca_cvSTE
#' @param bw Bandwidth for GWPCA.
#' @param x Data matrix.
#' @param loc Location matrix.
#' @param dMat Distance matrix.
#' @param k Number of neighbors for GWPCA.
#' @param robust Logical, indicating whether to use robust GWPCA.
#' @param kernel Kernel function for GWPCA.
#' @param adaptive Logical, indicating whether to use adaptive bandwidth.
#' @param cvContrib Logical, indicating whether to return a vector of
#' contributions for each location.
#' @param verbose Logical, indicating whether to display progress information.
#'
#' @returns If cvContrib is FALSE, returns the sum of contributions to the score
#' statistic used in cross-validation. If cvContrib is TRUE, returns a vector of
#' contributions for each location. The cvContrib == FALSE is used when a
#' bandwidth is automatically selected using `gwpca_bwSTE` function. The
#' cvContrib == TRUE is used when we need a score for each location to identify
#' discrepencies and outliers.
#'
#' @details The contribution of each observation to the score statistic used in
#' cross-validation for GWPCA is calculated. Outliers are taken to correspond to
#' high score (residual) values.
#'
#' @importFrom foreach %dopar%
#'
#' @seealso
#' \code{\link{wpca}}, \code{\link{rwpca}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Example usage
#' # Set parameters
#' bw <- 0.5
#' k <- 5
#'
#' # Create toy matrices
#' x <- matrix(rnorm(1000), ncol = 10)
#' loc <- abs(matrix(rnorm(200), ncol = 2))
#' dMat <- GWmodel::gw.dist(loc)
#'
#' # Calculate
#' gwpca_cvSTE(bw, x, loc, dMat, k)
#'
#' @export
gwpca_cvSTE <- function(bw,
                         x,
                         loc,
                         dMat,
                         k,
                         robust = FALSE,
                         kernel = "gaussian",
                         adaptive = FALSE,
                         cvContrib = FALSE,
                         verbose = FALSE) {

    ## Set handlers for a verbose approach
    if (verbose) {
        set.verbose(.global = TRUE,
                    .handlers = c("progress"))
    }

    ## Check for distance matrix
    if (missing(dMat)) {
        stop("A distance matrix needs to be provided with dMat argument!")
    }

    ## Select the type of WPCA
    if (robust == FALSE) {
        pcafun = wpca.ste
    } else {
        pcafun = rwpca.ste
    }

    ## Pre-allocate space
    n <- nrow(loc)
    score.contrib <- numeric(n)
    score <- 0

    ## Run CV contribution calculations
    if (verbose) {
        ## Run verbose
        score.contrib <- progrFUN(foreach = 1:n,
                                  FUN = int.gwpca.cv,
                                  ForEachArgs = list(int = 1:n, .combine = "c"),
                                  .bw = bw,
                                  .x = x,
                                  .dMat = dMat,
                                  .k = k,
                                  .kernel = kernel,
                                  .adaptive = adaptive,
                                  .pcafun = pcafun)

    } else {
        score.contrib <- foreach(int = 1:n, .combine = "c") %dopar% {
          int <- int # scope thing
          int.gwpca.cv(i = int,
                       .bw = bw,
                       .x = x,
                       .dMat = dMat,
                       .k = k,
                       .kernel = kernel,
                       .adaptive = adaptive,
                       .pcafun = pcafun)
        }
    }

    ## Return
    if (cvContrib) {
        return(score.contrib)
    } else {
        score <- sum(score.contrib)
        if (adaptive) {
          cat("Adaptive bandwidth(number of nearest neighbours):", bw,
              "CV score:", score, "\n")
        } else {
          cat("Fixed bandwidth:", bw, "CV score:", score, "\n")
        }

        return(score)
    }

}

# ---------------------------------------------------------------------------- #
#' Automatically find optimal bandwidth
#'
#' @name gwpca_bwSTE
#'
#' @description
#' A function to automatically find an optimal fixed or adaptive bandwidth.
#' This function is used to calibrate a basic or robust GWPCA via a
#' cross-validation approach.
#'
#' @param sfe A `SpatialFeatureExperiment` object.
#' @param vars A vector of variable names.
#' @param k The number of retained components. Default is 2.
#' @param robust If TRUE, apply robust GWPCA; otherwise, apply basic GWPCA.
#' Default is FALSE.
#' @param kernel The kernel function to use. Options include:
#'   - "gaussian": Gaussian kernel (default)
#'   - "exponential": Exponential kernel
#'   - "bisquare": Bisquare kernel
#'   - "tricube": Tricube kernel
#'   - "boxcar": Boxcar kernel
#' @param adaptive If TRUE, calculate an adaptive kernel where the bandwidth
#' corresponds to the number of nearest neighbors (adaptive distance). Default
#' is FALSE, where a fixed kernel is used (fixed distance).
#' @param p The order of the norm for Minkowski distance. If p = 1, it
#' represents Manhattan distance; if p = 2, it represents Euclidean distance.
#' Default is 2.
#' @param dMat The distance matrix. Default is NULL.
#'
#' @return Numeric, the optimal bandwidth.
#'
#' @importFrom SpatialFeatureExperiment spatialCoords
#' @importFrom SummarizedExperiment assay
#' @importFrom GWmodel gw.dist gold
#' @importFrom sp bbox
#' @importFrom magrittr %>%
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Load a SpatialFeatureExperiment object
#' data(sfe)
#'
#' # Prepare some parameters
#' vars <- rownames(sfe)[1:50]
#'
#' # Find the optimal bandwidth using basic GWPCA
#' \dontrun{
#' bw <- gwpca_bwSTE(sfe = sfe, vars = vars)
#' bw
#' }
#'
#' # Find the optimal bandwidth using robust GWPCA and an adaptive kernel
#' \dontrun{
#' bw <- gwpca_bwSTE(sfe = sfe, vars = vars,
#' robust = TRUE, adaptive = TRUE)
#' bw
#' }
#'
#' # Find the optimal bandwidth using a pre-calculated distance matrix
#' \dontrun{
#' dMat <- gw.dist(dp.locat = spatialCoords(sfe), p = 2)
#' bw <- gwpca_bwSTE(sfe = sfe, vars = vars, dMat = dMat)
#' bw
#' }
#'
#' @export
gwpca_bwSTE <- function(sfe,
                         vars,
                         k = 2,
                         robust = FALSE,
                         kernel = "gaussian",
                         adaptive = FALSE,
                         p = 2,
                         dMat) {

  ## Check sfe is an SFE object and extract some data
  if (is(sfe, "SpatialFeatureExperiment")) {
    dp.locat <- spatialCoords(sfe)
    data <- assay(sfe, assay) %>% t()

  } else {
    stop("Given data must be a SpatialFeatureExperiment object")
  }

  ## More checks
  if (missing(vars)) {
    stop("Variables input error: 'vars' argument must not be empty!")
  }

  ## Generate distance matrix
  dMat <- gw.dist(dp.locat = dp.locat, rp.locat = dp.locat, p = p,
                  theta = 0, longlat = FALSE)

  ## Extract the data that are going to be used
  col.nm <- colnames(data)
  var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if (length(var.idx) == 0) {
    stop("Variables input doesn't match with data")
  }

  x <- as.matrix(data[,var.idx]) # subset only for the required features
  # var.nms <- colnames(x)         # get feature names
  var.n <- ncol(x)               # get number of features
  dp.n <- nrow(dp.locat)         # get number of locations
  # ep.n <- nrow(elocat)           # get number of evaluation locations
  len.var <- length(vars)        # get length of 'vars' argument vector

  ## Check that all genes the user asked for are present in the selected data
  if (len.var > var.n) {
    stop("Invalid variables have been specified, please check them again!
         Not all variables provided in the 'vars' argument match the existing
         variables in the dataset provided.")
  }

  ## Find the range of the fixed bandwidth
  if (adaptive) {
    upper <- dp.n
    lower <- 2
  } else {
    if (!missing(dMat)) {
      upper <- range(dMat)[2]
      lower <- upper / 5000
    } else {
      dMat <- NULL
      if (p == 2) {
        b.box <- bbox(dp.locat)
        upper <-
          sqrt((b.box[1, 2] - b.box[1, 1])^2 + (b.box[2, 2] - b.box[2, 1])^2)
        lower <- upper / 5000
      } else {
        upper <- sapply(1:dp.n, function(i) {
          dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, p = p,
                             theta = 0, longlat = FALSE)
          range(dist.vi)[2]
        })
        upper <- max(upper)
        lower <- upper / 5000
      }
    }
  }

  ## Run the Golden selection optimisation algorithm
  bw <- NA
  bw <- gold(fun = gwpca_cvSTE,
             xL = lower,
             xU = upper,
             adapt.bw = adaptive,
             x = x,
             loc = dp.locat,
             k = k,
             robust = robust,
             kernel = kernel,
             adaptive = adaptive,
             p = p,
             dMat = dMat)

  return(bw)

}

# -----------------------------------------------------------------------------#
#' Print GWPCA output
#'
#' @name print.gwpca
#' @description
#' A function to print GWPCA results in the prompt
#'
#' @param x An object of class 'gwpca'.
#' @param ... Not used currently. Left for future additions.
#'
#' @importFrom stats printCoefmat
#' @importFrom methods as
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @export
print.gwpca <- function(x, ...) {
  if (!inherits(x, "gwpca")) {
    stop("It's not a gwpca object")
  }

  cat("   ******************************************************************\n")
  cat("   *                       Package   GWmodel                        *\n")
  cat("   ******************************************************************\n")
  cat("   Call:\n")
  cat("   ")
  vars <- x$GW.arguments$vars
  cat("\n   Variables concerned: ", vars)
  cat("\n   The number of retained components: ", x$GW.arguments$k)
  dp.n <- dim(x$loadings)[1]
  cat("\n   Number of data points:", dp.n)

  cat("\n   ****************************************************************\n")
  cat("   *                Results of Principal Components Analysis        *\n")
  cat("   ******************************************************************\n")
  print(summary(x$pca, loadings = TRUE, cutoff = 0))

  cat("\n   ****************************************************************\n")
  cat("   *    Results of Geographically Weighted Prin. Comp. Analysis     *\n")
  cat("   ******************************************************************\n")
  cat("\n   *******************Model calibration information****************\n")
  cat("   Kernel function for geographically weighting:", x$GW.arguments$kernel,
      "\n")

  if (x$GW.arguments$adaptive) {
    cat("   Adaptive bandwidth for geographically and temporally weighting: ",
        x$GW.arguments$bw, " (number of nearest neighbours)\n", sep = "")
  } else {
    cat("   Fixed bandwidth for geographically and temporally weighting: ",
        x$GW.arguments$bw, "\n")
  }

  if (x$GW.arguments$p == 2) {
      cat("   Distance metric for geographically weighting: Euclidean \n")
  } else if (x$GW.arguments$p == 1) {
      cat("   Distance metric for geographically weighting: Manhattan \n")
  } else if (is.infinite(x$GW.arguments$p)) {
      cat("   Distance metric for geographically weighting: Chebyshev \n")
  } else {
      cat("   Distance metric for geographically weighting: A generalized
          Minkowski distance metric is used with p=", x$GW.arguments$p, ".\n")
  }

  cat("\n   ************     Summary of GWPCA information:    **************\n")
  var.names <- paste("Comp", 1:x$GW.arguments$k, sep = ".")
  cat("   Local variance: \n")
  local.SD <- t(apply(x$var[, 1:x$GW.arguments$k], 2, summary))[, c(1:3, 5, 6)]
  rownames(local.SD) <- var.names
  printCoefmat(local.SD)

  cat("   Local Proportion of Variance: \n")
  local.PV <- t(apply(as(x$SDF, "data.frame")[, 1:(x$GW.arguments$k + 1),
                                              drop = FALSE],
                      2,
                      summary))[, c(1:3, 5, 6)]
  rownames(local.PV) <- c(var.names, "Cumulative")
  printCoefmat(local.PV)

  cat("\n   ****************************************************************\n")

  invisible(x)
}

# ----------------------------NON-EXPORTED FUNCTIONS-------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Progress function for foreach loop.
#'
#' @name progrFUN
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#'
#' @keywords internal
#'
#' @param foreach Numerical vector to iterate over.
#' @param FUN Function to apply to each iteration.
#' @param ForEachArgs Additional arguments to pass to the foreach loop.
#' @param ... Additional arguments passed to the FUN function.
#'
#' @importFrom progressr progressor
#'
#' @details This function is used as a progress function for the foreach loop.
#' It displays the progress of the loop iteration using the progressor function.
#'
#' @return Returns the result of applying the FUN function to each iteration.
#'
#' @seealso
#' \code{\link{progressor}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
progrFUN <- function(foreach, FUN, ForEachArgs = list(), ...) {
  ## Set progress handler
  pr <- progressr::progressor(along = foreach)

  ## Run the foreach
  do.call("foreach", c(ForEachArgs)) %dopar% {
    ## Fetch progress
    pr()
    ## Run the function
    int <- int # scope thing
    FUN(i = int, ...)
  }
}
# ---------------------------------------------------------------------------- #
#' A basic WPCA with distance weighting.
#'
#' @name wpca.ste
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' Performs a basic Weighted-PCA with distance weighting.
#'
#' @param x Data matrix.
#' @param wt Weight vector.
#' @param ... Additional arguments passed to the svd function.
#'
#' @keywords internal
#'
#' @return Returns the result of the singular value decomposition.
#'
#' @details
#' This function performs a basic Weighted Principal Component Analysis (WPCA)
#' with distance weighting. The gene counts are centred by subtracting a
#' distance-weighted mean, and then the centred matrix is multiplied by the
#' square root of the weights. Finally, the Singular Value Decomposition (SVD)
#' is applied to the weighted and centred matrix.
#'
#' @seealso
#' \code{\link{svd}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
wpca.ste <- function(x, wt, ...) {
  ## Centre the gene counts by subtracting a distance-weighted mean
  x_centered <- t(t(x) - base::colSums(x * wt) / sum(wt))
  ## Multiply the centred matrix by the square root of the weights
  x_weighted <- x_centered * sqrt(wt)
  ## Perform Single Value Decomposition
  svd(x_weighted, ...)
}
# ---------------------------------------------------------------------------- #
#' A robust SVD function.
#'
#' @name robustSvd.ste
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' Performs a robust Singular Value Decomposition (SVD) using the MCD estimator.
#'
#' @param x Data matrix.
#' @param alpha Tuning parameter for the MCD estimator (default = 0.75).
#'
#' @keywords internal
#'
#' @return Returns a list with the loadings (v) and variances (d) obtained from
#' the SVD.
#'
#' @details
#' This function performs a robust Singular Value Decomposition (SVD) using the
#' Minimum Covariance Determinant (MCD) estimator. The MCD estimator is used to
#' estimate the covariance matrix of the data, and then the SVD is applied to
#' the estimated covariance matrix. The loadings (v) and variances (d) obtained
#' from the SVD are returned. The true variances (d) are returned, unlike in the
#' basic GWPCA case, i.e., d=(pc$sdev)^2.
#'
#' @importFrom stats princomp
#' @importFrom robustbase covMcd
#'
#' @seealso
#' \code{\link{princomp}}, \code{\link{covMcd}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
robustSvd.ste <- function(x, alpha = 0.75) {
  cov_matrix <- covMcd(x, alpha = alpha)$cov
  pc <- princomp(covmat = cov_matrix)
  return(list(v = pc$loadings, d = pc$sdev))
}
# ---------------------------------------------------------------------------- #
#' Calculate a distance weighted median for use in Robust WPCA.
#'
#' @name wt.median.ste
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' Calculates a weighted median for use in Robust Weighted Principal Component
#' Analysis (WPCA) with distance weighting.
#'
#' @keywords internal
#'
#' @param x Data matrix.
#' @param wt Weight vector.
#'
#' @return Returns a vector of weighted medians, one for each column of the
#' input matrix.
#'
#' @details
#' This function calculates a weighted median for each column of the input data
#' matrix. The weighted median is computed by rearranging the data in ascending
#' order, computing the cumulative sum of weights, and identifying the index of
#' the median. The weighted median is then returned for each column.
#'
#' @seealso
#' \code{\link{order}}
#'
## (note using medians as a robust estimate of location)
### Calculate medians
wt.median.ste <- function(x, wt) {
  wt.median.1.ste <- function(.x, .wt) {
    ox <- order(.x)                    # rearrange in ascending order.
    wox <- cumsum(.wt[ox])             # compute the cumulative sum of weights.
    posn <- which.min(abs(wox - 0.5))  # identifies the index of the median.
    return(.x[ox][posn])
  }
  apply(x, 2, wt.median.1.ste, .wt = wt)
}
# ---------------------------------------------------------------------------- #
#' A robust WPCA.
#'
#' @name rwpca.ste
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' Performs a robust Weighted Principal Component Analysis (WPCA) with distance
#' weighting.
#'
#' @keywords internal
#'
#' @param x Data matrix.
#' @param wt Weight vector.
#' @param nu Number of observations to use as initial estimates (default = 0).
#' @param nv Number of components to retain (default = 2).
#'
#' @return Returns a list with the loadings (v) and variances (d) obtained from
#' the WPCA.
#'
#' @details
#' This function performs a robust Weighted Principal Component Analysis (WPCA)
#' with distance weighting. It calculates the medians of each column of the data
#' matrix and subtracts them from the data to centre it. The centred matrix is
#' then multiplied by the weight vector to apply the distance weighting. The
#' resulting weighted matrix is subjected to a robust Singular Value
#' Decomposition (SVD) using the `robustSvd.ste` function. The loadings (v) and
#' variances (d) obtained from the SVD are returned, with the number of
#' components retained determined by the `nv` parameter.
#'
#' @importFrom robustbase covMcd
#' @importFrom stats princomp
#'
#' @seealso
#' \code{\link{robustSvd.ste}}, \code{\link{wt.median.ste}},
#' \code{\link{princomp}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
### Weighted PCA
rwpca.ste <- function(x, wt, nu = 0, nv = 2) {
  medians <- wt.median.ste(x, wt)
  centered <- t(t(x) - medians)
  weighted_centered <- centered * wt
  result <- robustSvd.ste(weighted_centered)
  result$v <- result$v[, 1:nv]
  return(result)
}
# ---------------------------------------------------------------------------- #
#' GWPCA calculating for loop
#'
#' @name int.gwpca.FLoop
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' Internal function used in the Geographically Weighted Principal Component
#' Analysis (GWPCA).
#'
#' @keywords internal
#'
#' @param i The index of the current location.
#' @param ..pcafun The PCA function to use.
#' @param ..x The data matrix.
#' @param ..dMat The distance matrix.
#' @param ..k The number of components to retain.
#' @param ..bw The bandwidth parameter.
#' @param ..adaptive Logical indicating if adaptive bandwidth should be used.
#' @param ..kernel The kernel function to use for weighting.
#' @param ..scores Logical indicating if local scores should be calculated.
#'
#' @details
#' This internal function calculates the GWPCA for a single location (index i)
#' using the specified PCA function. It selects the weights for the location
#' based on the distance matrix and bandwidth parameter. The function skips the
#' calculation if the bandwidth is too small (less than or equal to 5). If the
#' calculation proceeds, it runs the PCA function on the subset of data with
#' non-zero weights, and stores the resulting weights and scores (if specified).
#'
#' @returns The GWPCA outcome as a list.
#'
#' @importFrom GWmodel gw.weight
#'
#' @seealso
#' \code{\link{gwpca}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
int.gwpca.FLoop <- function(i,
                            ..pcafun,
                            ..x,
                            ..dMat,
                            ..k,
                            ..bw,
                            ..adaptive,
                            ..kernel,
                            ..scores) {
  ## Select weights for location i
  wt <- gw.weight(vdist = ..dMat[,i], bw = ..bw,
                  kernel = ..kernel, adaptive = ..adaptive)
  use <- wt > 0

  ## Skip if the bandwidth is too small
  if (length(wt) <= 5) {
    warning(paste("Too small bandwidth at location: ", i, " and the results
      can't be given there."))
  }

  ## Run GWPCA
  temp <- ..pcafun(..x[use,], wt[use], nu = 0, nv = ..k)

  ## Store weights
  temp$wt <- wt

  ## Calculate the local scores using matrix operations
  if (..scores) {
    temp$score <- ..x[use,] %*% temp$v
  } else {
    temp$score <- NULL
  }

  ## Return GWPCA outcome
  return(temp)
}
# ---------------------------------------------------------------------------- #
#' GWPCA core
#'
#' @name int.gwpca
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' Internal function used in the Geographically Weighted Principal Component
#' Analysis (GWPCA).
#'
#' @keywords internal
#'
#' @param .sfe The spatial feature experiment object.
#' @param .assay The name of the assay to use from the object.
#' @param .elocat The evaluation locations.
#' @param .vars The variables to include in the analysis.
#' @param .p The power parameter for distance calculation.
#' @param .k The number of components to retain.
#' @param .bw The bandwidth parameter.
#' @param .kernel The kernel function to use for weighting.
#' @param .adaptive Logical indicating if adaptive bandwidth should be used.
#' @param .scores Logical indicating if local scores should be calculated.
#' @param .robust Logical indicating if robust WPCA should be used.
#' @param .cv Logical indicating if cross-validation should be performed.
#' @param .verbose Logical indicating if verbose output should be displayed.
#'
#' @details
#' This internal function performs the GWPCA analysis using the specified
#' parameters. It takes a spatial feature experiment object (.sfe) and extracts
#' the necessary data. The evaluation locations (.elocat) are checked, and if
#' not provided or in an invalid format, the default locations from the spatial
#' feature experiment object are used. The function also performs checks on the
#' variables, bandwidth, and other input parameters. It generates the distance
#' matrix and selects the subset of data for analysis based on the variables.
#' The global PCA is run using the selected data. The memory is pre-allocated
#' for storing the results of the loop. The type of WPCA (robust or non-robust)
#' is determined based on the input parameter. The GWPCA is then performed
#' either in parallel or in a verbose mode. The results are stored in the
#' appropriate arrays and lists. Cross-validation is performed if requested, and
#' the results are added to the output. The influence of weights or variances in
#' the data is adjusted if not using robust WPCA. The final output includes the
#' PCA results, loadings, spatial points data frame (SDF), GWPCA scores
#' (if calculated), variance, local percent of variation, input arguments,
#' cross-validation results, and geometry information.
#'
#' @returns A list of class "gwpca".
#'
#' @importFrom stats princomp
#' @importFrom foreach foreach
#' @importFrom Matrix rowSums
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom GWmodel gw.dist
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialFeatureExperiment spatialCoords
#' @importFrom DelayedArray DelayedArray
#'
#' @seealso
#' \code{\link{gwpcaSTE}}, and the original function from GWmodel package;
#' [gwpca()].
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
int.gwpca <- function(.sfe,
                      .assay,
                      .elocat,
                      .vars,
                      .p,
                      .k,
                      .bw,
                      .kernel,
                      .adaptive,
                      .scores,
                      .robust,
                      .cv,
                      .verbose) {

  ## Check sfe is an SFE object and extract some data
  if (is(.sfe, "SpatialFeatureExperiment")) {
    dp.locat <- spatialCoords(.sfe)

    data <- SummarizedExperiment::assay(.sfe, .assay)
    ## Check matrix class
    if (is(data, "dgCMatrix")) {
      data <- DelayedArray::DelayedArray(data) %>% t()
    } else {
      data <- t(data)
    }

  } else {
    stop("Given data must be a SpatialFeatureExperiment object")
  }

  ## Check for evaluation locations
  if (is.null(.elocat)) {
    .elocat <- dp.locat

  } else if (is.numeric(.elocat) && dim(.elocat)[2] == 2) {
    .elocat <- .elocat

  } else {
    warning("Output loactions are not a two-column numeric vector")
    .elocat <- dp.locat
  }

  ## More checks
  if (missing(.vars)) {
    stop("Variables input error: 'vars' argument must not be empty!")
  }

  if (missing(.bw) || .bw <= 0) {
    stop("Bandwidth is specified incorrectly: either missing or is <= 0")
  }

  ## Generate distance matrix
  dMat <- gw.dist(dp.locat = dp.locat, rp.locat = .elocat, p = .p,
                  theta = 0, longlat = FALSE)

  ## Extract the data that are going to be used
  col.nm <- colnames(data)
  var.idx <- match(.vars, col.nm)[!is.na(match(.vars, col.nm))]
  if (length(var.idx) == 0) {
    stop("Variables input doesn't match with data")
  }

  x <- as.matrix(data[,var.idx]) # subset only for the required features
  var.nms <- colnames(x)         # get feature names
  var.n <- ncol(x)               # get number of features
  dp.n <- nrow(dp.locat)         # get number of locations
  ep.n <- nrow(.elocat)          # get number of evaluation locations
  len.var <- length(.vars)       # get length of 'vars' argument vector

  ## Check that all genes the user asked for are present in the selected data
  if (len.var > var.n) {
    stop("Invalid variables have been specified, please check them again!
         Not all variables provided in the 'vars' argument match the existing
         variables in the dataset provided.")
  }

  ## Run a global PCA
  pca.res <- princomp(x, cor = TRUE, scores = .scores)

  ## Pre-allocate memory for the loop
  w <- array(data = NA, c(ep.n, var.n, .k))
  d <- array(data = NA, c(ep.n, var.n))

  gwpca.scores <- NULL
  if (.scores) {
    gwpca.scores <- vector("list", ep.n)
  }

  ## Select the type of WPCA
  if (.robust == FALSE) {
    pcafun = wpca.ste
  } else {
    pcafun = rwpca.ste
  }

  ## Run GWPCA ----------------------------------- ##
  if (.verbose) {
    ## Run verbose
    temp <- progrFUN(foreach = 1:ep.n,
                     FUN = int.gwpca.FLoop,
                     ForEachArgs = list(int = 1:ep.n),
                     ..pcafun = pcafun,
                     ..x = x,
                     ..dMat = dMat,
                     ..k = .k,
                     ..bw = .bw,
                     ..adaptive = .adaptive,
                     ..kernel = .kernel,
                     ..scores = .scores)

  } else {
    temp <- foreach::foreach(int = 1:ep.n) %dopar% {
      int <- int # scope thing
      int.gwpca.FLoop(i = int,
                      ..pcafun = pcafun,
                      ..x = x,
                      ..dMat = dMat,
                      ..k = .k,
                      ..bw = .bw,
                      ..adaptive = .adaptive,
                      ..kernel = .kernel,
                      ..scores = .scores)
    }
  }

  for (ep in 1:ep.n) {
    ## Assign values to w
    w[ep,,] <- temp[[ep]]$v
    ## Assign values to d
    d[ep,] <- temp[[ep]]$d
    ## Assign values to gwpca.scores
    if (.scores) {
      gwpca.scores[[ep]] <- temp[[ep]]$score
    }
    ## Assign values to wt
    wt <- temp[[ep]]$wt
  }

  ## Add dimension names to array
  if (!is.null(rownames(x))) {
    ## Row names (location names)
    dimnames(w)[[1]] <- rownames(x)
  }

  if (!is.null(colnames(x))) {
    ## Column names (feature names)
    dimnames(w)[[2]] <- colnames(x)
  }

  dimnames(w)[[3]] <- paste0("PC", 1:.k) # PC names

  ## --------------------------------------- ##
  ## Perform Cross-Validation (CV) for selected bandwidth (bw)
  ## Pre-set a vector of length equal to number of locations (dp.n)
  CV <- list()

  ## Check if Cross-Validation is asked and run it
  if (.cv) {
    message("Performing Cross-Validation for selected bandwith ", .bw)
    CV$CV <- gwpca_cvSTE(x = x,
                         loc = dp.locat,
                         dMat = dMat,
                         bw = .bw,
                         k = .k,
                         robust = .robust,
                         kernel = .kernel,
                         adaptive = .adaptive,
                         cvContrib = TRUE,
                         verbose = .verbose)

    outliers <- outlierCutoff(CV$CV, coef = 1.5)
    CV$out_up <- outliers$out_up
    CV$out_down <- outliers$out_down
    CV$is_outlier <- outliers$outliers

    ## Add arguments to a list for the output
    GW.arguments <- list(vars = .vars,
                         k = .k,
                         bw = .bw,
                         kernel = .kernel,
                         adaptive = .adaptive,
                         p = .p,
                         dp.n = dp.n,
                         scores = .scores)

  }

  ## ------------------------------------------- ##
  ## Adjust for the influence of weights or variances in the data if not robust
  d1 <- if (.robust) d^2 else (d / (sum(wt)^0.5))^2
  local.PV <- d1[, 1:.k] / rowSums(d1) * 100     # get local percent of variation
  var.names <- paste("Comp", 1:.k, sep = ".")    # get the Prin. Comp. names
  win.var.pc1 <- max.col(abs(w[,,1]))   # get the gene with max loading per spot
  res.df <- data.frame(local.PV, rowSums(local.PV), .vars[win.var.pc1])
  names(res.df) <- c(paste(var.names, "PV", sep = "_"),
                     "local_CP",
                     "win_var_PC1")
  SDF <- SpatialPointsDataFrame(coords = .elocat,
                                data = res.df,
                                match.ID = FALSE)
  geometry <- colGeometry(.sfe, "spotHex")
  ## Put them in a list to output
  res <- list(pca = pca.res, loadings = w, SDF = SDF,
              gwpca.scores = gwpca.scores, var = d1,
              local.PV = local.PV, GW.arguments = GW.arguments,
              CV = CV, geometry = geometry)
  class(res) <- "gwpca"

  invisible(res)

}
# ---------------------------------------------------------------------------- #
#' Internal: Corss Validation
#'
#' @name int.gwpca.cv
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' Internal function used in the Geographically Weighted Principal Component
#' Analysis (GWPCA) cross-validation.
#'
#' @keywords internal
#'
#' @param i The index of the evaluation location.
#' @param .bw The bandwidth parameter.
#' @param .x The data matrix.
#' @param .dMat The distance matrix.
#' @param .k The number of components to retain.
#' @param .kernel The kernel function to use for weighting.
#' @param .adaptive Logical indicating if adaptive bandwidth should be used.
#' @param .pcafun The PCA function to use.
#'
#' @details
#' This internal function performs cross-validation for a given evaluation
#' location.
#' It takes the index of the evaluation location (i) and calculates the weights
#' using the specified bandwidth (.bw), distance matrix (.dMat), kernel function
#' (.kernel), and adaptive bandwidth flag (.adaptive). The weight for the
#' evaluation location (i) is set to 0 to exclude it from the calculation.
#' The function checks if the number of neighboring locations with positive
#' weights is greater than 1, and if not, returns an infinite value for the
#' cross-validation result and displays a warning message. The PCA function
#' (.pcafun) is applied to the subset of data with positive weights,
#' and the eigenvectors are obtained. The eigenvectors are then used to
#' calculate the cross-validation result based on the squared difference
#' between the original data at the evaluation location (i) and the data
#' reconstructed using the eigenvectors.
#'
#' @returns The cross-validation result as a vector.
#'
#' @importFrom GWmodel gw.weight
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
int.gwpca.cv <- function(i,
                         .bw,
                         .x,
                         .dMat,
                         .k,
                         .kernel,
                         .adaptive,
                         .pcafun) {
  wt <- gw.weight(.dMat[,i], .bw, .kernel, .adaptive)
  wt[i] <- 0
  use <- wt > 0

  if (sum(use) <= 1) {
    out <- Inf
    out <- Inf
    cat("Too small bandwidth:", .bw, "and the CV value can't be given there.\n")
    return(NULL)
  }

  v <- .pcafun(.x[use, ], wt[use], nu = 0, nv = .k)$v
  v <- v %*% t(v)

  out <- sum((.x[i, ] - .x[i, ] %*% v))^2
}
