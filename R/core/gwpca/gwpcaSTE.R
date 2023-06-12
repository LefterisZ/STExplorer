#' @name gwpca.ste
#' 
#' @description
#' A Geographically Weighted Principal Components Analysis function based on 
#' \code(gwpca) function from \code{GWmodel} package. The function is re-written
#' and optimised to work faster and with a \code{SpatialFeatureExperimetn} (SFE)
#' object.
#' 
#' 
#' 

#' --------------------------------------------------------------------------- #
#' @name porgrFUN
#' @param foreach numerical vector to iterate over
#'
progrFUN <- function(foreach, FUN, ForEachArgs = list(), ...) {
  ## Set progress handler
  pr <- progressor(along = foreach)
  
  ## Run the foreach
  do.call("foreach", c(ForEachArgs)) %dopar% {
    ## Fetch progress
    pr()
    ## Run the function
    FUN(i = i, ...)
  }
}
#' --------------------------------------------------------------------------- #
#' @name wpca
#' @description
#' A basic WPCA with distance weighting.
#' 
wpca.ste <- function(x, wt, ...) {
  ## Centre the gene counts by subtracting a distance-weighted mean
  x_centered <- t(t(x) - colSums(x * wt) / sum(wt))
  ## Multiply the centred matrix by the square root of the weights
  x_weighted <- x_centered * sqrt(wt)
  ## Perform Single Value Decomposition
  svd(x_weighted, ...)
}
#' --------------------------------------------------------------------------- #
#' @name robustSvd
#' @description
#' A robust SVD function.
#' 
## The true variances (d) are returned, unlike in the basic GWPCA case 
## i.e., d=(pc$sdev)^2
robustSvd.ste <- function(x, alpha = 0.75) {
  cov_matrix <- covMcd(x, alpha = alpha)$cov
  pc <- princomp(covmat = cov_matrix)
  return(list(v = pc$loadings, d = pc$sdev))
}
#' --------------------------------------------------------------------------- #
#' @name wt.median
#' @description
#' A function to calculate a weighted median for use in Robust WPCA with 
#' distance weighting.
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
  apply(x, 2, wt.median.1, .wt = wt)
}
#' --------------------------------------------------------------------------- #
#' @name rwpca
#' @description
#' A robust WPCA with distance weighting function.
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
#' --------------------------------------------------------------------------- #
#' @name int.gwpca.FLoop
#' 
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
    next
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
#' --------------------------------------------------------------------------- #
#' @name int.gwpca
#' 
#' 

int.gwpca <- function(.obj, 
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
  
  ## Check obj is an SFE object and extract some data
  if (is(.obj, "SpatialFeatureExperiment")) {
    dp.locat <- spatialCoords(.obj)
    data <- assay(.obj, .assay) %>% t()
    
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
                     ForEachArgs = list(i = 1:ep.n),
                     ..pcafun = pcafun,
                     ..x = x,
                     ..dMat = dMat,
                     ..k = .k,
                     ..bw = .bw,
                     ..adaptive = .adaptive,
                     ..kernel = .kernel,
                     ..scores = .scores)
    
  } else { 
    temp <- foreach::foreach(i = 1:ep.n) %dopar% { 
      int.gwpca.FLoop(i = i,
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
  
  for (i in 1:ep.n) {
    ## Assign values to w
    w[i,,] <- temp[[i]]$v
    ## Assign values to d
    d[i,] <- temp[[i]]$d
    ## Assign values to gwpca.scores
    if (.scores) {
      gwpca.scores[[i]] <- temp[[i]]$score
    }
    ## Assign values to wt
    wt <- temp[[i]]$wt
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
    
    outliers <- outlier_cutoff(CV$CV, coef = 1.5)
    CV$out_up <- outliers$out_up
    CV$out_down <- outliers$out_down
    CV$is_outlier <- outliers$outliers
    
    ## Add arguments to a list for the output
    GW.arguments <- list(vars = vars,
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
  geometry <- colGeometry(.obj, "spotHex")
  ## Put them in a list to output
  res <- list(pca = pca.res, loadings = w, SDF = SDF, 
              gwpca.scores = gwpca.scores, var = d1, 
              local.PV = local.PV, GW.arguments = GW.arguments, 
              CV = CV, geometry = geometry)
  class(res) <- "gwpca"
  
  invisible(res)
  
}
# ---------------------------------------------------------------------------- #
#' @name int.gwpca.cv
#' 
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
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' @name gwpca_ste
#' 
#' @description
#' A Geographically Weighted Principal Components Analysis function based on 
#' \code(gwpca) function from \code{GWmodel} package. The function is re-written
#' and optimised to work faster and with a \code{SpatialFeatureExperimetn} (SFE)
#' object.
#' 
#' @param obj a \code{SpatialFeatureExperiment} (SFE) object.
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
#' gaussian: wgt = exp(-.5*(vdist/bw)^2);
#' exponential: wgt = exp(-vdist/bw);
#' bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
#' tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#' boxcar: wgt=1 if dist < bw, wgt=0 otherwise
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
#' @export

gwpcaSTE <- function(obj, 
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
  locations <- nrow(spatialCoords(obj))
  message("Locations to iterate over: ", locations)
  message("Please be patient. This will take a while...")
  mins <- round((locations * 2) / 60, 3)
  message("E.T.A. is approximatelly: ", mins, " minutes")
  message("Sit back and relax!! :)")
  ## fetch start time
  s <- Sys.time()
  
  ## RUN GWPCA
  gwpca_result <- int.gwpca(.obj = obj, 
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
  print(gwpca_result)
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
#' @name gwpca.cv.ste
#' @description
#' @param cvContrib defaults to FALSE and gives back a single score. When used
#' inside \code{gwpca} function is set to TRUE and gives back a vector of scores
#' one for each location.
#' 
# Contribution of each observation to the score statistic used in cross-
# validation for gwpca
# Outliers taken to correspond to high score (residual) values... 
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
        pcafun = wpca
    } else {
        pcafun = rwpca
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
                                  ForEachArgs = list(i = 1:n, .combine = "c"),
                                  .bw = bw,
                                  .x = x,
                                  .dMat = dMat,
                                  .k = k,
                                  .kernel = kernel,
                                  .adaptive = adaptive,
                                  .pcafun = pcafun)
        
    } else {
        score.contrib <- foreach(i = 1:n, .combine = "c") %dopar% {
          int.gwpca.cv(i = i,
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
#' @name bw.gwpca.ste
#' @description
#' A Function to automatically find an optimal fixed or adaptive bandwidth. A 
#' function for automatic bandwidth selection to calibrate a basic or robust 
#' GWPCA via a cross-validation approach only.
#' 
#' @export

gwpca_bwSTE <- function(obj,
                         vars,
                         k = 2,
                         robust = FALSE,
                         kernel = "gaussian",
                         adaptive = FALSE,
                         p = 2,
                         dMat) {
  
  ## Check obj is an SFE object and extract some data
  if (is(obj, "SpatialFeatureExperiment")) {
    dp.locat <- spatialCoords(obj)
    data <- assay(obj, assay) %>% t()
    
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
  var.nms <- colnames(x)         # get feature names
  var.n <- ncol(x)               # get number of features
  dp.n <- nrow(dp.locat)         # get number of locations
  ep.n <- nrow(elocat)           # get number of evaluation locations
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
  bw <- gold(fun = gwpca.cv.ste,
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

#' ----------------------------------------------------------------------------#
#' @name print.gwpca
#' @description
#' A function to print GWPCA results in the prompt
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

