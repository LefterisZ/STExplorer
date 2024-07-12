# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Load packages ----
library(pryr)
library(lineprof)
library(microbenchmark)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Prepare input data ----

# data
# pop
# distmat
# kind = NA
# ncluster = 2
# m = 2
# distance = "euclidean"
# order = 2
# alpha = 0.7
# a = 1
# b = 1
# max.iter = 500
# error = 1e-05
# randomN = 0
# uij = NA
# vi = NA

# fgwc_in
# k_range = 2:10
# index_type = "composite"
# elbow_method = "knee"
# m_sfe
# sample_id = NULL
# algorithm = "classic"
# distMat = NULL
# dMetric = NULL
# parameters = fgwc_params(algorithm = "classic", ncluster = 5)
# plot = TRUE
# verbose = FALSE

# m_sfe
sfe_test <- sfe[, 1:1000]
sfe_nmf_test <- sfe_nmf[1:1000,]
sfe_test <- addDistMat(sfe_test, 2)
sample_id = NULL
data = sfe_nmf_test
pop = rep(1, nrow(data))
distMat = getDistMat(sfe_test, dMetric = "euclidean", sample_id = NULL)
dMetric = NULL
algorithm = "classic"
parameters = fgwc_params(algorithm = "classic", ncluster = 5)

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Re-implemented FGWC source code ----
## Check SFE or MSFE?
sfe1 <- .int_sfeORmsfe(m_sfe = sfe_test, sample_id = sample_id)
c(address(sfe1), refs(sfe1))
c(address(sfe_test), refs(sfe_test))

## Get distance Matrix
if (is.null(distMat)) {
  distMat <- getDistMat(sfe1, dMetric = dMetric, sample_id = sample_id)
}

## Get population
if (is.null(pop)) {
  pop <- rep(1, nrow(data))
}

## Set defaults
if (is.null(parameters) && algorithm == "classic") {
  parameters <- c(kind = 'u', ncluster = 2, m = 2,
                  distance = 'euclidean', order = 2,
                  alpha = 0.7, a = 1, b = 1,
                  max.iter = 500, error = 1e-5, randomN = 1)
}

main_params <- list(data = data,
                    pop = pop,
                    distmat = getDistMat(sfe_test, dMetric = "euclidean", sample_id = NULL))

## Run FGWC
if (algorithm == "classic") {
  ## Run classic FGWC
  fgwc <- do.call(naspaclust::fgwcuv, c(parameters, main_params))
} else if (algorithm == "abc") {
  ## Run FGWC with Artificial Bee Colony (ABC) optimisation
  fgwc <- do.call(naspaclust::abcfgwc, c(parameters, main_params))
}

## Add extra information in the fgwc output
fgwc$finaldata <- .int_addToFGWCout(fgwc = fgwc,
                                    sfe = sfe,
                                    geom_type = "hex")

fgwc$metageneSignatures <- attr(data, "basis")


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


data = data
pop = pop
distmat = getDistMat(msfe = sfe_test, dMetric = "euclidean")
kind = "u"
ncluster = 5
m = 1.5
distance = "manhattan"
order = 1
alpha = 0.5
a = 1.8
b = 1
max.iter = 500
error = 1e-05
randomN = 1
uij = NA
vi = NA

ptm <- proc.time()
stopifnot(kind == "v" || kind == "u" || is.na(kind) == TRUE)

## Ensure 'data' is a matrix
if (!is.matrix(data)) {
  data <- as.matrix(data)
}

## Set up required variables
n <- nrow(data)
d <- ncol(data)
beta <- 1 - alpha

iter <- 0
conv <- numeric(0)

## Handle alpha = 1 case
if (alpha == 1) {
  pop <- rep(1, n)
  distmat <- matrix(1, n, n)
}

## Ensure 'pop' is a matrix with one column
pop <- matrix(pop, ncol = 1)
mi.mj <- pop %*% t(pop)

## Initialise 'uij' membership matrix
if (is.na(kind) || kind == "u") {
  if (is.na(uij)) { # If 'uij' is NA, we need to initialise it.
    ## Set the random seed to ensure reproducibility.
    set.seed(randomN)
    ## Generate a uniform random matrix of dimensions n x ncluster.
    uij <- matrix(runif(n * ncluster), n, ncluster)
    ## Normalise each row of 'uij' so that the sum of each row equals 1.
    ## This makes 'uij' a proper fuzzy membership matrix.
    uij <- uij / rowSums(uij)
  } else {
    ## If 'uij' is not NA, it means it's externally provided.
    ## Convert it to a matrix of dimensions n x ncluster.
    uij <- matrix(uij, n, ncluster)
  }
}

## Main loop for FGWC "u" (membership approach)
## Commentary PER LINE:
##  - Initialise a matrix 'old_uij' with the same dimensions as 'uij'
##    (n x ncluster), filled with zeros. It is used to store the membership
##    matrix from the previous iteration to monitor convergence.
##  - Start the iterative process.
##    Continue until the maximum change in 'uij' is less than 'error' or the
##    maximum number of iterations 'max.iter' is reached.
##  - Update 'old_uij' to keep the previous state of 'uij' to compare changes
##    between iterations.
##  - Compute the new centroids 'vi' using the current state of 'uij'.
##    New centroids vi are a weighted mean of the data points, where the
##    weights are given by the membership degrees old_uij^m. 't(old_uij^m)' is
##    the transposed matrix of 'old_uij' raised element-wise to the power of
##    'm'. This matrix is multiplied by 'data', and then element-wise divided
##    by the column sums of 'old_uij^m' to get the new centroids
##  - Compute the distances from each data point to each centroid,
##    adjust them by the fuzzifier 'm', and then take element-wise inversion.
##    'rdist::cdist(data, vi, distance, order)' computes the distance between
##    each pair of data points and centroids. After squaring and adjusting by
##    'm', the distances are used to update 'uij'.
##  - Update 'uij' by taking the inverse of the row sums of the inverse of
##    'dists'. This ensures that each element of 'uij' is between 0 and 1 and
##    all elements in a row sum to 1.
##  - Modify 'uij' using the spatial and population information.
##    --> This step incorporates the geographical component into the clustering.
##  - Increment the iteration counter.
##  - Calculate the convergence metric for this iteration and append it to the
##    'conv' vector. This uses the updated 'uij' and the centroids 'vi' to
##    compute the fuzzy within-cluster sum of squares.
if (is.na(kind) || kind == "u") { # FGWC "u" (membership approach)
  ## Initialise a matrix 'old_uij' with the same dimensions as 'uij'
  old_uij <- matrix(0, n, ncluster)
  ## Start the iterative process.
  while (max(abs(uij - old_uij)) > error && iter < max.iter) {
    ## Update 'old_uij' to keep the previous state of 'uij' to compare changes between iterations.
    old_uij <- uij
    ## Compute the new centroids 'vi' using the current state of 'uij'.
    vi <- (t(old_uij^m) %*% data) / colSums(old_uij^m)
    ## Compute adjusted by fuzzifier 'm' distances from each data point to each centroid
    dists <- (rdist::cdist(data, vi, distance, order)^2)^(1 / (m - 1))
    ## Update 'uij'. Ensure range between 0-1 and a row sum to 1 for each element of 'uij'.
    uij <- (1 / dists) / rowSums(1 / dists)
    ## Modify 'uij' using the spatial and population information.
    uij <- .int_membershipUpdate(data, uij, mi.mj, distmat, alpha, beta, a, b)
    ## Increment the iteration counter.
    iter <- iter + 1
    ## Calculate the convergence metric.
    conv <- c(conv, sum(uij^m * (rdist::cdist(data, vi, distance, order)^2)))
  }
}

## Main loop for FGWC "v" (centroid approach)
## Commentary PER LINE:
##  - Check if 'vi' (the centroid matrix) is not initialised (is NA).
##    If true, initialise it.
##  - Generate initial centroids 'vi' using the 'gen_vi' function with
##    uniform distribution.
##  - Initialise 'v_new' with the current centroids 'vi'.
##  - Initialise 'v_old' as an array of 'Inf' with the same dimensions as 'vi'
##    to ensure the loop starts.
##  - Initialise the membership matrix 'uij' with zeros, with 'n' rows and
##    'ncluster' columns.
##  - Start the iterative process; continue until centroids change less than
##    'error' or maximum iterations 'max.iter' are reached.
##  - Update 'v_old' to keep the previous centroids to compare changes between
##    iterations.
##  - Compute the distances from each data point to each centroid 'v_old',
##    adjust by the fuzzifier 'm', and take element-wise inversion.
##    'rdist::cdist(data, v_old, distance, order)' computes the distance
##    between each pair of data points and centroids. After squaring and
##    adjusting by 'm', the distances are used to update 'uij'.
##  - Update 'uij' by taking the inverse of the row sums of the inverse of
##    'dists'. This ensures that each element of 'uij' is between 0 and 1 and
##    all elements in a row sum to 1.
##  - Modify 'uij' using the spatial and population information.
##    This step incorporates the geographical component into the clustering.
##  - Compute the new centroids 'v_new' using the current state of 'uij'.
##    't(uij^m)' is the transposed matrix of 'uij' raised element-wise to the
##    power of 'm'. This matrix is multiplied by 'data', and then element-wise
##    divided by the column sums of 'uij^m' to get the new centroids.
##  - Calculate the convergence metric for this iteration and append it to the
##    'conv' vector. This uses the updated 'uij' and the centroids 'v_new' to
##    compute the fuzzy within-cluster sum of squares.
##  - Increment the iteration counter.
##  - Update the final centroids 'vi' with 'v_new' after finishing the iterations.
if (kind == "v") { # FGWC "v" (centroid approach)
  if (is.na(vi)) {
    ## Generate initial centroids.
    vi <- .int_generateInitialCentroids(data, ncluster, "uniform", randomN)
  }
  ## Initialise 'v_new' with the current centroids 'vi'.
  v_new <- vi
  ## Initialise 'v_old' as an array of 'Inf' with the same dimensions as 'vi'.
  v_old <- array(Inf, dim = dim(vi))
  ## Initialise the membership matrix 'uij' with zeros.
  dists <- (rdist::cdist(data, v_new, distance, order)^2)^(1 / (m - 1))
  uij <- (1 / dists) / rowSums(1 / dists)
  ## Start the iterative process.
  while (max(abs(v_new - v_old)) > error && iter < max.iter) {
    ## Update 'v_old' to keep the previous centroids to compare changes between iterations.
    v_old <- v_new
    ## Compute adjusted by the fuzzifier 'm' distances from each data point to each centroid 'v_old'.
    dists <- (rdist::cdist(data, v_old, distance, order)^2)^(1 / (m - 1))
    ## Update 'uij'. Ensure range between 0-1 and a row sum to 1 for each element of 'uij'.
    uij <- (1 / dists) / rowSums(1 / dists)
    ## Modify 'uij' using the spatial and population information.
    uij <- .int_membershipUpdate(data, uij, mi.mj, distmat, alpha, beta, a, b)
    ## Compute the new centroids 'v_new' using the current state of 'uij'.
    v_new <- (t(uij^m) %*% data) / colSums(uij^m)
    ## Calculate the convergence metric for this iteration.
    conv <- c(conv, sum(uij^m * (rdist::cdist(data, v_new, distance, order)^2)))
    ## Increment the iteration counter.
    iter <- iter + 1
  }
  ## Update the final centroids 'vi' with 'v_new' after finishing the iterations.
  vi <- v_new
}

fgwc_obj <- sum(uij^m * (rdist::cdist(data, vi, distance, order)^2))
finaldata <- .int_determineCluster(data, uij)
cluster <- finaldata[, ncol(finaldata)]

## Prepare output
result <- list(
  converg = conv,
  f_obj = fgwc_obj,
  membership = uij,
  centroid = vi,
  validation = .int_getIndexes(data, cluster, uij, vi, m, exp(1)),
  iteration = iter,
  cluster = cluster,
  finaldata = finaldata,
  call = match.call(),
  time = proc.time() - ptm
)
class(result) <- 'fgwc'

## Clean result object
result$call$data <- NULL
result$call$pop <- NULL
result$call$distmat <- NULL

## Output
printOut <- c(order, ncluster, m, randomN)
names(printOut) <- c("order", "ncluster", "m", "randomN")
print(printOut)
return(result)

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Re-implemented Index functions ----
# Set the seed for reproducibility
set.seed(42)

# Generate a matrix of data points (20 data points, 2 features)
data <- matrix(rnorm(40), nrow = 20, ncol = 2)

# Generate a fuzzy membership matrix uij (20 data points, 3 clusters)
uij <- matrix(runif(40), nrow = 20, ncol = 2)
# Normalize each row to sum to 1
uij <- uij / rowSums(uij)

# Generate centroids vi (3 centroids, 2 features)
vi <- matrix(rnorm(4), nrow = 2, ncol = 2)

# Fuzzifier m
m <- 2

# Clusters
fg <- naspaclust::fgwcuv(data, rep(1, 20), distmat = distmat[1:20,1:20])
cluster <- fg$finaldata$cluster

# Show the generated data
print("Data (first 5 rows):")
print(head(data, 5))

print("Membership uij (first 5 rows):")
print(head(uij, 5))

print("Centroids vi:")
print(vi)

## Separation Index Original
SI1 <- function(data, uij, vi, m) {
  d <- matrix(0, nrow(data), nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i, j] <- sum((vi[j,] - data[i,])^2)
    }
  }
  vkvi <- matrix(0, nrow(vi), nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i, k] <- sum((vi[k,] - vi[i,])^2)
    }
  }
  diag(vkvi) <- Inf
  pt1 <- sum((uij^2) * d)
  pt2 <- nrow(data) * min(vkvi)
  return(sum(pt1 / pt2))
}

## Separation Index optimised
optimisedSI1 <- function(data, uij, vi, m) {
  d <- as.matrix(dist(rbind(data, vi)))
  d <- d[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(vi))]^2

  viDist <- as.matrix(dist(vi))^2
  diag(viDist) <- Inf

  pt1 <- sum((uij^2) * d)
  pt2 <- nrow(data) * min(viDist)

  return(pt1 / pt2)
}

## Original XB1 Function
XB1 <- function(data, uij, vi, m) {
  d <- matrix(0, nrow(data), nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i, j] <- sum((vi[j,] - data[i,])^2)
      if (d[i, j] == 0) {
        d[i, j] <- Inf
      }
    }
  }
  pt1 <- sum((uij^m) * d)
  pt2 <- nrow(data) * min(d)
  return(pt1 / pt2)
}

## Optimised XB1 Function
optimisedXB1 <- function(data, uij, vi, m) {
  d <- as.matrix(dist(rbind(data, vi)))
  d <- d[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(vi))]^2
  d[d == 0] <- Inf

  pt1 <- sum((uij^m) * d)
  pt2 <- nrow(data) * min(d)

  return(pt1 / pt2)
}

## Original SC1 function
SC1 <- function(data,cluster,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  pt1 <- colSums((uij^m)*d)
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[i,]-vi[k,])^2)
    }
    Ni <- length(which(cluster==i))
    vkvi[i,]<-Ni*vkvi[i,]
  }
  pt2 <- colSums(vkvi)
  return(sum(pt1/pt2))
}

## Optimised SC1 function
optimisedSC1 <- function(data, cluster, uij, vi, m) {
  d <- as.matrix(dist(rbind(data, vi)))
  d <- d[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(vi))]^2
  pt1 <- colSums((uij^m) * d)
  viDist <- as.matrix(dist(vi))^2

  # Step 2: Precompute the number of points in each cluster
  # Ni will have a length of max(cluster) to handle all possible cluster indices
  Ni <- numeric(max(cluster))
  for (i in 1:max(cluster)) {
    Ni[i] <- sum(cluster == i)
  }

  # Step 3: Initialize vkvi and scale rows by the corresponding Ni
  vkvi <- viDist
  for (i in 1:nrow(vi)) {
    if (i <= length(Ni)) {
      vkvi[i,] <- Ni[i] * vkvi[i,]
    }
  }

  pt2 <- colSums(vkvi)

  return(sum(pt1 / pt2))
}

## Original IFV1 function
IFV1 <- function(data,uij,vi,m) {
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[k,]-vi[i,])^2)
    }
  }
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  sigmaD <- sum(d)/(nrow(data)*nrow(vi))
  SDmax <- max(vkvi)
  log2u <- colSums(log(uij,2))/nrow(data)
  u2ij <- colSums(uij^2)
  inside <- sum(u2ij*(log(nrow(vi),2)-log2u))
  return(sum(u2ij*((log(nrow(vi),2)-log2u)^2)/nrow(data)*(SDmax/sigmaD)))
}

## Optimised IFV1 function
optimisedIFV1 <- function(data, uij, vi, m) {
  ## Calculate vkvi matrix using vectorized operations
  distance_matrix <- as.matrix(dist(vi, method = "euclidean"))^2
  vkvi <- distance_matrix

  ## Calculate d matrix using vectorized operations
  d <- as.matrix(dist(rbind(data, vi), method = "euclidean"))^2
  d <- d[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(vi))]

  ## Calculate sigmaD
  sigmaD <- sum(d) / (nrow(data) * nrow(vi))

  ## Calculate SDmax
  SDmax <- max(vkvi)

  ## Calculate log2u and u2ij
  log2u <- colSums(log(uij, 2)) / nrow(data)
  u2ij <- colSums(uij^2)

  ## Compute inside term
  log_term <- log(nrow(vi), 2) - log2u
  inside <- sum(u2ij * log_term)

  ## Calculate final IFV value
  result <- sum(u2ij * ((log_term)^2) / nrow(data) * (SDmax / sigmaD))

  return(result)
}

## Original Kwon1 function
Kwon1 <- function(data,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  s <- matrix(0,nrow(vi))
  vivj <- matrix(0,nrow(vi),nrow(vi))
  for (j in 1:nrow(vi)) {
    s[j,] <- (sum(vi[j,]-colMeans(data))^2)
  }
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  for (i in 1:nrow(vi)) {
    for (j in 1:nrow(vi)) {
      vivj[i,j] <- sum((vi[i,]-vi[j,])^2)
    }
  }
  diag(vivj) <- Inf
  pt1 <- colSums((uij^m)*d)
  pt2 <- sum(s)/nrow(vi)
  pt3 <- min(vivj)
  return(sum((pt1+pt2)/pt3))
}

## Optimised Kwon1 function
optimisedKwon1 <- function(data, uij, vi, m) {
  ## Calculate d matrix using vectorized operations
  d <- as.matrix(dist(rbind(data, vi), method = "euclidean"))^2
  d <- d[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(vi))]

  ## Calculate s vector using vectorized operations
  s <- colSums((vi - colMeans(data))^2)

  ## Calculate vivj matrix using vectorized operations
  vivj <- as.matrix(dist(vi, method = "euclidean"))^2
  diag(vivj) <- Inf  # Set diagonal elements to Inf

  ## Calculate pt1, pt2, and pt3
  pt1 <- colSums((uij^m) * d)
  pt2 <- sum(s) / nrow(vi)
  pt3 <- min(vivj)

  ## Compute the final Kwon value
  result <- sum((pt1 + pt2) / pt3)

  return(result)
}

## Test the original and optimised SI1 functions
SI1_result <- SI1(data, uij, vi, m)
optimisedSI1_result <- optimisedSI1(data, uij, vi, m)

cat("SI1 Result:", SI1_result, "\n")
cat("Optimised SI1 Result:", optimisedSI1_result, "\n")

## Test the original and optimised XB1 functions
XB1_result <- XB1(data, uij, vi, m)
optimisedXB1_result <- optimisedXB1(data, uij, vi, m)

cat("XB1 Result:", XB1_result, "\n")
cat("Optimised XB1 Result:", optimisedXB1_result, "\n")

## Test the original and optimised SC1 functions
SC1_result <- SC1(data, uij, vi, m)
optimisedSC1_result <- optimisedSC1(data, uij, vi, m)

cat("SC1 Result:", SC1_result, "\n")
cat("Optimised SC1 Result:", optimisedSC1_result, "\n")

## Test the original and optimised IFV1 functions
IFV1_result <- IFV1(data, uij, vi, m)
optimisedIFV1_result <- optimisedIFV1(data, uij, vi, m)

cat("IFV1 Result:", IFV1_result, "\n")
cat("Optimised IFV1 Result:", optimisedIFV1_result, "\n")

## Test the original and optimised Kwon1 functions
Kwon1_result <- Kwon1(data, uij, vi, m)
optimisedKwon1_result <- optimisedKwon1(data, uij, vi, m)

cat("Kwon1 Result:", Kwon1_result, "\n")
cat("Optimised Kwon1 Result:", optimisedKwon1_result, "\n")

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Test running times ----
mbm <- microbenchmark("STE" = fgwcSTE(m_sfe = sfe_test,
                                      data = data,
                                      pop = pop,
                                      distMat = distMat,
                                      parameters = parameters),
                      "original" = naspaclust::fgwc(data = data,
                                                    pop = pop,
                                                    distmat = distMat,
                                                    fgwc_param = parameters))


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Testing C++ running times ----
ptm1 <- proc.time()
res <- t(uij)
proc.time() - ptm1

ptm2 <- proc.time()
res <- Rfast::transpose(uij)
proc.time() - ptm2

# --------------------------------------- #
ptm1 <- proc.time()
res <- rowSums(uij)
proc.time() - ptm1

ptm2 <- proc.time()
res <- Rfast::rowsums(uij)
proc.time() - ptm2

# --------------------------------------- #
uij <- matrix(runif(n * ncluster), n, ncluster)
uij <- uij / rowSums(uij)
old_uij <- uij
m = 1.5
data = sfe_nmf

ptm1 <- proc.time()
res <- (t(old_uij^m) %*% data) / colSums(old_uij^m)
proc.time() - ptm1

ptm2 <- proc.time()
res2 <- (Rfast::transpose(old_uij^m) %*% data) / Rfast::colsums(old_uij^m)
proc.time() - ptm2

# --------------------------------------- #
vi = res
distance = "euclidean"
roder = 2

ptm1 <- proc.time()
res <- (rdist::cdist(data, vi, distance, order)^2)^(1 / (m - 1))
proc.time() - ptm1

ptm2 <- proc.time()
res2 <- (Rfast::transpose(old_uij^m) %*% data) / Rfast::colsums(old_uij^m)
proc.time() - ptm2

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Housekeeping ----
rm(k_range, index_type, elbow_method, sample_id, algorithm, distMat, dMetric,
   parameters, plot, verbose, data, pop, main_params, distMat, fgwc, kind,
   ncluster, m, distance, order, alpha, a, b, max.iter, error, randomN, uij,
   vi, n, d, beta, conv, iter, .int_generateInitialCentroids,
   .int_membershipUpdate, .int_determineCluster, dists, mi.mj, new_uij, old_uij,
   uij, v_new, v_old, vi, vi2, vi3, vi4, sfe1, maxs, means, mins, sds, SI1,
   optimisedSI1, XB1, optimisedXB1, optimisedSC1, SC1, optimisedIFV1, IFV1,
   Kwon1, optimisedKwon1, Kwon1_result, optimisedIFV1_result, XB1_result,
   optimisedKwon1_result, optimisedSI1_result, optimisedXB1_result, IFV1_result,
   SI1_result)
