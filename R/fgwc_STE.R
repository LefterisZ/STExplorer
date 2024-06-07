#' Identify the optimal number of NMF factors
#'
#' This function determines the optimal number of factors for Non-negative
#' Matrix Factorization (NMF) by evaluating the reconstruction error of the
#' input matrix `A` and a randomized version of `A`. The evaluation is done
#' using either Mean Squared Error (MSE) or Mean Kullback-Leibler (MKL)
#' divergence.
#'
#' @param m_sfe The \code{SpatialFeatureExperiment} or
#' \code{MetaSpatialFeatureExperiment} object containing spatial expression
#' data.
#' @param sample_id The sample ID or index for which to perform clustering. If
#' \code{NULL}, all samples are used.
#' @param assay The name of the assay in \code{m_sfe} to use for clustering
#' (default is "logcounts").
#' @param top_hvgs A character vector of gene names representing the top highly
#' variable genes. The gene names are expected to be ENSGene IDs.
#' @param k_range A numeric vector specifying the range of factor numbers to
#' evaluate (default is \code{seq(2, 12, 2)}).
#' @param n_cores An integer indicating the number of cores to use for parallel
#' processing (default is 1).
#' @param do_plot Logical, indicating whether to plot the reconstruction error
#' results (default is \code{TRUE}).
#' @param seed An integer specifying the random seed for reproducibility
#' (default is 1).
#' @param loss Character, either \code{"mse"} (mean squared error) or
#' \code{"mkl"} (mean Kullback-Leibler divergence) (default is \code{"mse"}).
#' @param max.iter An integer specifying the maximum number of iterations for
#' NMF (default is 250).
#'
#' @return An integer. The number of optimal factors.
#'
#' @details
#' This function is based on Frigyesi et al, 2008.
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2623306/ and on the \code{swne}
#' package by Yan Wu https://doi.org/10.1016/j.cels.2018.10.015
#' (https://github.com/yanwu2014/swne/).
#' This function evaluates the reconstruction error for a specified range of
#' factors (\code{k_range}) to determine the optimal number of factors for NMF.
#' It compares the reconstruction error of the input matrix and a randomized
#' version to identify the optimal cut-off.
#'
#' @examples
#' \dontrun{
#' # Assuming `sfe` is your SpatialFeatureExperiment object
#' optimal_factors <- fgwc_nmfFactorNumber(m_sfe = sfe,
#'                                         sample_id = "ID1",
#'                                         top_hvgs = top_genes)
#' }
#'
#' @seealso
#' \code{\link[scater]{calculateNMF}} for additional details on NMF clustering.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords clustering spatial-expression nmf
#' @family clustering functions
#' @rdname fgwc_nmfFactorNumber
#' @aliases fgwc_nmfFactorNumber
#'
#' @export
fgwc_nmfFactorNumber <- function(m_sfe,
                                 sample_id = NULL,
                                 assay = "logcounts",
                                 top_hvgs,
                                 k_range = seq(2, 10, 1),
                                 n_cores = 1,
                                 do_plot = TRUE,
                                 seed = 1,
                                 loss = c("mse", "mkl"),
                                 max.iter = 250) {
  ## SFE or metaSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get data
  nmf_input <- assay(sfe, assay)[rownames(sfe) %in% top_hvgs,]
  message("Initialising...")

  ## Perform Checks
  .int_initialise(A = nmf_input,
                  k_range = k_range,
                  loss = loss,
                  seed = seed)

  A <- as.matrix(nmf_input)
  A.rand <- matrix(sample(A), nrow(A), ncol(A))
  message("Finding optimal number...")

  ## Find optimum number of factors
  result <- .int_findOptimalFactors(A = A,
                                    A.rand = A.rand,
                                    k_range = k_range,
                                    n_cores = n_cores,
                                    max.iter = max.iter,
                                    loss = loss,
                                    seed = seed)

  ## Print if required
  if (do_plot) {
    print(plotFGWC_factorSelection(result, font.size = 12))
  }

  return(result)
}


#' Perform Non-negative Matrix Factorisation (NMF)
#'
#' This function performs Non-negative Matrix Factorisation (NMF) on the gene
#' expression matrix present inside a SpatialFeatureExperiment or a
#' MetaSpatialFeatureExperiment.
#'
#' @param m_sfe The \code{SpatialFeatureExperiment} or
#' \code{MetaSpatialFeatureExperiment} object containing spatial expression
#' data.
#' @param sample_id The sample ID or index for which to perform clustering.
#' @param assay The name of the assay in \code{m_sfe} to use for clustering
#' (default is "logcounts").
#' @param top_hvgs A character vector of gene names representing the top highly
#' variable genes. The gene names are expected to be ENSGene IDs.
#' @param ncomponents The number of components (clusters) to identify
#' (default is 2). The components are the 'factors' calculated mby NMF.
#' @param ntop The maximum number of top features (genes) to use in the
#' clustering (default is 600).
#' @param subset_row An optional vector specifying the subset of rows to use in
#' the clustering.
#' @param scale Logical, indicating whether to scale the input data
#' (default is FALSE).
#' @param ... Additional arguments to be passed to \code{scater::calculateNMF}.
#'
#' @return An NMF result object with factors representing the reduced dimensions
#'  of the original expression matrix.
#'
#' @examples
#' \dontrun{
#' # Assuming `sfe` is your SpatialFeatureExperiment object
#' nmf_result <- fgwc_nmf(m_sfe = sfe, sample_id = "ID1", top_hvgs = top_genes)
#' }
#'
#' @seealso
#' \code{\link[scater]{calculateNMF}} for additional details on NMF clustering.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords clustering spatial-expression nmf
#' @family clustering functions
#' @rdname fgwc_nmf
#' @aliases fgwc_nmf
#' @importFrom dplyr select
#'
#' @export
fgwc_nmf <- function(m_sfe,
                     sample_id,
                     assay = "logcounts",
                     top_hvgs,
                     ncomponents,
                     ntop = 600,
                     subset_row = NULL,
                     scale = FALSE,
                     ...) {
  ## SFE or metaSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get data
  nmf_input <- assay(sfe, assay)[rownames(sfe) %in% top_hvgs,]

  ## Run NMF
  nmf <- .int_calculate_nmf(x = nmf_input,
                            ncomponents = ncomponents,
                            ntop = ntop,
                            subset_row = subset_row,
                            scale = scale,
                            seed = 1,
                            ...)
  colnames(nmf) <- sprintf("Factor%02d", 1:ncol(nmf))
  rownames(nmf) <- colnames(sfe)

  return(nmf)
}


#' Generate Parameters for Fuzzy Geographically Weighted Clustering
#'
#' This function generates parameters for Fuzzy Geographically Weighted
#' Clustering (FGWC) or optimization algorithms for FGWC.
#'
#' @param algorithm The type of FGWC or optimization algorithm to use. Choose
#' from "classic", "abc", "fpa", "gsa", "hho", "ifa", "pso", "tlbo".
#' @param ... arguments passed down to the specific algorithm function
#'
#' @details
#' depending on the algorithm the function calls internally the respective
#' function.
#'
#' @return A list of FGWC parameters or optimization algorithm
#' parameters.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' params <- fgwc_params(algorithm = "classic", ncluster = 3, kind = "u", m = 2)
#' }
#'
#' @seealso \code{\link[naspaclust]{fgwcuv}}, \code{\link[naspaclust]{abcfgwc}},
#' \code{\link[naspaclust]{fpafgwc}}, \code{\link[naspaclust]{gsafgwc}},
#' \code{\link[naspaclust]{hhofgwc}}, \code{\link[naspaclust]{ifafgwc}},
#' \code{\link[naspaclust]{psofgwc}}, \code{\link[naspaclust]{tlbofgwc}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords clustering fuzzy-logic optimization fgwc
#' @family clustering functions
#' @rdname fgwc_params
#' @aliases fgwc_params
#'
#' @export
fgwc_params <- function(algorithm = c("classic", "abc", "fpa", "gsa",
                                      "hho", "ifa", "pso", "tlbo"),
                        ...) {
  ## Check algorithm arguments
  algorithm <- match.arg(algorithm)

  if (algorithm == "classic") {
    ## A named vector that consists of FGWC parameter, or a vector that consists
    ##   of optimization algorithm parameter.
    params <- classic_params(...)
    return(params)
  } else if (algorithm == "abc") {
    ## optimization algorithm parameter.
    params <- abc_params(...)
    return(params)
  }
  else if (!algorithm %in% c("classic", "abc")) {
    stop("\nCurrently this function supports only the 'classic' algortihm",
         "type and the 'abc' optimisation method.\n",
         "For the other optimisation algorithms, refer to the vignette with\n",
         "the default parameters for any optimisation algortihm or refer to\n",
         "the help page (?fgwc_params) to find the optimisation algorithm's\n",
         "help under section 'See also'.")
  }

}

#' @param ncluster The number of clusters.
#' @param kind Use \code{'u'} for the membership approach and \code{'v'} for
#' the centroid approach.
#' @param m Degree of fuzziness or fuzzifier (default is 2).
#' @param distance The distance metric between data and centroid (default is
#' "euclidean").
#' @param order Minkowski order (default is 1).
#' @param alpha The old membership effect with `[0,1]`. If \code{alpha} equals
#' 1, it will be the same as Fuzzy C-Means. If 0, it equals the neighborhood
#' effect.
#' @param a Spatial magnitude of distance (default is 1).
#' @param b Spatial magnitude of population (default is 1).
#' @param max.iter Maximum iteration (default is 500).
#' @param error Error tolerance (default is 1e-5).
#' @param randomN Random seed for initialization (if \code{uij} or \code{vi}
#' is NA, default is 1).
#' @param uij Membership matrix initialization.
#' @param vi Centroid matrix initialization.
#'
#' @rdname fgwc_params
#' @export
classic_params <- function(ncluster,
                           kind = "u",
                           m = 1.5,
                           distance = "manhattan",
                           order = 1,
                           alpha = 0.5,
                           a = 2,
                           b = 1,
                           max.iter = 500,
                           error = 1e-5,
                           randomN = 1,
                           uij = NA,
                           vi = NA) {
  ## A named vector that consists of FGWC parameter
  params <- list(kind = kind,
                 ncluster = ncluster,
                 m = m,
                 distance = distance,
                 order = order,
                 alpha = alpha,
                 a = a, # Modifying the membership matrix utilizing the distance matrix.
                 b = b, # Modifying the membership matrix utilizing the population.
                 max.iter = max.iter,
                 error = error,
                 randomN = randomN,
                 uij = uij,
                 vi = vi)
  return(params)
}

#' @param error error tolerance. Default is 1e-5.
#' @param max.iter maximum iteration. Default is 500.
#' @param randomN random seed for initialisation (if uij or vi is NA).
#' Default is 0.
#' @param vi.dist a string of centroid population distribution between
#' "uniform" (default) and "normal". Can be defined as vi.dist= in opt_param.
#' @param nfood number of foods population. Can be defined as npar= in
#' opt_param.
#' @param n.onlooker number of onlooker bees, Can be defined as n.onlooker in
#' opt_param.
#' @param limit number of turns to eliminate food with no solutions. Can be
#' defined as limit in opt_param.
#' @param pso whether to add PSO term in bee's movement. Either TRUE or FALSE.
#' Can be defined as pso in opt_param.
#' @param abc.same number of consecutive unchanges to stop the iteration. Can
#' be defined as same= in opt_param. description
#' @rdname fgwc_params
#' @export
abc_params <- function(ncluster,
                       m = 1.6,
                       distance = "manhattan",
                       order = 1,
                       alpha = 0.5,
                       a = 2,
                       b = 1,
                       error = 1e-05,
                       max.iter = 100,
                       randomN = 0,
                       vi.dist = "uniform",
                       nfood = 10,
                       n.onlooker = 5,
                       limit = 4,
                       pso = FALSE,
                       abc.same = 10) {
  ## optimization algorithm parameter.
  params <- list(ncluster = ncluster,
                 m = m,
                 distance = distance,
                 order = order,
                 alpha = alpha,
                 a = a, # Modifying the membership matrix utilizing the distance matrix.
                 b = b, # Modifying the membership matrix utilizing the population.
                 error = error,
                 max.iter = max.iter,
                 randomN = randomN,
                 vi.dist = vi.dist,
                 nfood = nfood,
                 n.onlooker = n.onlooker,
                 limit = limit,
                 pso = pso,
                 abc.same = abc.same)
  return(params)
}

#' Identify the Optimal Number of Clusters Using FGWC
#'
#' This function identifies the optimal number of clusters (`k`) for Fuzzy
#' Geographically Weighted Clustering (FGWC) by evaluating different clustering
#' indices over a range of `k` values and finding the elbow point.
#'
#' @param fgwc_in The matrix (locations x factors/components) that will be used
#' for the actual FGWC as generated by the \code{fgwc_nmf} function.
#' @param k_range An integer vector specifying the range of `k` values to
#' evaluate (default is \code{2:10}).
#' @param index_type A character string indicating which index type to use for
#' evaluation. Must be one of \code{"composite"}, \code{"FPC"}, \code{"CE"},
#' \code{"SC"}, or \code{"XB"} (default is \code{"FPC"}).
#' @param elbow_method A character string specifying the method for detecting
#' the elbow point. Options are \code{"angle"} or \code{"knee"} (default is
#' \code{"angle"}).
#' @param m_sfe The \code{SpatialFeatureExperiment} or
#' \code{MetaSpatialFeatureExperiment} object containing spatial expression
#' data.
#' @param sample_id A character string indicating the name of the sample. It is
#' required when the `m_sfe` argument is provided with a MetaSFE object.
#' @param dMetric Character string specifying the distance metric.
#' @param distMat An `n*n` distance matrix between regions.
#' @param algorithm A character string indicating the algorithm used for FGWC.
#' Currently accepting one of these two options: \code{"classic"} or
#' \code{"abc"}.
#' @param parameters A vector that consists of FGWC parameters
#' (see \code{\link{fgwc_params}} for parameter details).
#' @param plot Logical, indicating whether to plot the index values with the
#' optimal `k` highlighted (default is \code{TRUE}).
#' @param verbose Logical, indicating whether to print progress messages
#' (default is \code{FALSE}).
#' @param ... Additional arguments passed to the \code{fgwcSTE} function.
#'
#' @return An integer representing the optimal number of clusters.
#'
#' @details
#' The function evaluates different clustering indices over a range of `k`
#' values and identifies the elbow point using either the angle or knee method.
#'
#' If the `index_type` is set to \code{"composite"}, a composite score is
#' calculated based on multiple indices, and the elbow point is determined
#' using the composite score.
#'
#' @examples
#' \dontrun{
#' # Assuming `fgwc_input` is a matrix generated by the `fgwc_nmf` function
#' # and `sfe` is your SpatialFeatureExperiment object
#' optimal_k <- fgwc_findOptimumK(fgwc_in = fgwc_input, m_sfe = sfe)
#' }
#'
#' @seealso \code{\link{fgwcSTE}}, \code{\link{fgwc_params}}
#' @keywords clustering spatial-expression fgwc
#' @rdname fgwc_findOptimumK
#' @export
fgwc_findOptimumK <- function(fgwc_in,
                              k_range = 2:10,
                              index_type = "FPC",
                              elbow_method = "angle",
                              m_sfe,
                              sample_id = NULL,
                              algorithm = "classic",
                              distMat = NULL,
                              dMetric = NULL,
                              parameters,
                              plot = TRUE,
                              verbose = FALSE,
                              ...) {
  ## Validate the index_type argument
  index_type <- toupper(index_type) # make sure is all caps
  if (!index_type %in% c("COMPOSITE", "FPC", "CE", "SC", "XB")) {
    stop("index_type must be one of 'composite', 'FPC', 'CE', 'SC', or 'XB'")
  }

  ## SFE or metaSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get distance Matrix for FGWC
  if (is.null(distMat)) {
    distMat <- getDistMat(sfe, dMetric = dMetric, sample_id = sample_id)
  }

  ## Initialise a data frame to store k and the corresponding index values
  if (index_type == "COMPOSITE") {
    index_data <- data.frame(k = integer(0),
                             FPC = integer(0),
                             SC = integer(0),
                             CE = integer(0),
                             XB = integer(0))
  } else {
    index_data <- data.frame(k = integer(0),
                             IndexValue = numeric(0),
                             IndexType = character(0))
  }

  ## Iterate over all k values in k_range
  for (k in k_range) {
    if (verbose) {
      message("Running FGWC for k = ", k, "...\n")
    }

    ## Update number of clusters in parameters list
    parameters$ncluster <- k

    ## Perform fuzzy geographically weighted clustering
    fgwc_result <- fgwcSTE(m_sfe = sfe,
                           sample_id = sample_id,
                           data = fgwc_in,
                           distMat = distMat,
                           dMetric = dMetric,
                           algorithm = algorithm,
                           parameters = parameters,
                           ...)

    ## Extract the validation indices
    indices <- fgwc_result$validation

    ## Determine the current index value based on the selected index_type
    if (index_type == "CE") {
      current_index_value <- indices$CE
    } else if (index_type == "SC") {
      current_index_value <- indices$SC
    } else if (index_type == "XB") {
      current_index_value <- indices$XB
    } else if (index_type == "FPC") {
      current_index_value <- .int_fpcoeff(fgwc_result$membership)
    }

    ## Append the current index value to the index_data data frame
    if (index_type == "COMPOSITE") {
      verbose <- FALSE
      current_index_value <- .int_fpcoeff(fgwc_result$membership)
      index_data <- rbind(index_data, data.frame(k = k,
                                                 FPC = current_index_value,
                                                 SC = indices$SC,
                                                 CE = indices$CE,
                                                 XB = indices$XB))
    } else {
      index_data <- rbind(index_data, data.frame(k = k,
                                                 IndexValue = current_index_value,
                                                 IndexType = index_type))
    }

    if (verbose) {
      message("k = ", k, ": ", index_type, " = ", current_index_value, "\n")
    }
  }

  ## If composite score is selected, calculate the score and update the index df
  if (index_type == "COMPOSITE") {
    index_data <- .int_compositeScore(index_data = index_data)
  }

  ## Find best k by calculating the elbow point
  if (elbow_method == "angle") {
    best_k <- .int_findElbowPointAngle(index_data = index_data)
    subtitle = "Largest Angle Between Pairs of k"
  } else if (elbow_method == "knee") {
    best_k <- .int_findElbowPointKnee(index_data = index_data)
    subtitle = "Knee (Largest Dist From Straight Line)"
  }

  ## Plot the index values if plot is TRUE
  if (plot) {
    p <- ggplot(data = index_data, aes(x = k, y = IndexValue)) +
      geom_line() +
      geom_point() +
      geom_point(data = index_data[index_data$k == best_k, ],
                 aes(x = k, y = IndexValue),
                 shape = 21, colour = "black", fill = "#69b3a2", size = 6) +
      geom_vline(xintercept = best_k, linetype = "dashed", color = "red") +
      labs(title = paste0("Optimal Clustering Evaluation with ", index_type, " Index"),
           subtitle = paste0("Elbow Method: ", subtitle),
           x = "Number of Clusters k",
           y = paste0(index_type, " Index Value")) +
      theme_classic() +
      theme(axis.text.x = element_text(hjust = 1, size = 15, colour = "black"),
            axis.text = element_text(size = 15, colour = "black"),
            axis.title = element_text(size = 15, colour = "black"),
            legend.position = "none")

    print(p)
  }

  ## Delete after
  # assign(paste0("index_data", index_type), index_data, envir = .GlobalEnv)

  ## Return the best k
  return(best_k)
}

#' Fuzzy Geographicaly Weighted Clustering
#'
#' @description Wrapper around classic fuzzy clustering with addition of
#' spatial configuration of membership matrix. It supports optimisation
#' algorithms.
#'
#' @param m_sfe The \code{SpatialFeatureExperiment} or
#' \code{MetaSpatialFeatureExperiment} object containing spatial expression
#' data.
#' @param sample_id A character string indicating the name of the sample. It is
#' required when the `m_sfe` argument is provided with a MetaSFE object.
#' @param data an object of data with d>1. Can be \code{matrix} or
#' \code{data.frame}. If your data is univariate, bind it with \code{1} to get
#' a data frame with 2 columns. We suggest this to be the dimensionality
#' reduction matrix (locations x factors/components) as generated by the
#' \code{fgwc_nmf} function.
#' @param pop an n*1 vector contains population.
#' @param dMetric Character string specifying the distance metric.
#' @param distMat an n*n distance matrix between regions.
#' @param algorithm algorithm used for FGWC. Currently accepting one of these
#' two options: \code{"classic"} or \code{"abc"}
#' @param parameters a vector that consists of FGWC parameters
#' (see \code{\link{fgwc_params}} for parameter details)
#'
#' @return an object of class \code{"fgwc"}.\cr
#' An \code{"fgwc"} object contains as follows:
#' \itemize{
#' \item \code{converg} - the process convergence of objective function
#' \item \code{f_obj} - objective function value
#' \item \code{membership} - membership matrix
#' \item \code{centroid} - centroid matrix
#' \item \code{validation} - validation indices (there are partition
#' coefficient (\code{PC}), classification entropy (\code{CE}),
#' SC index (\code{SC}), separation index (\code{SI}), Xie and Beni's index
#' (\code{XB}), IFV index (\code{IFV}), and Kwon index (Kwon))
#' \item \code{max.iter} - Maximum iteration
#' \item \code{cluster} - the cluster of the data
#' \item \code{finaldata} - The final data (with the cluster)
#' \item \code{call} - the syntax called previously
#' \item \code{time} - computational time.
#' }
#'
#' @details Fuzzy Geographically Weighted Clustering (FGWC) was developed by
#' \code{naspaclust} by adding neighbourhood effects and population to
#' configure the membership matrix in Fuzzy C-Means. There are two kinds of
#' options in doing classical FGWC. The first is using \code{"u"} (default)
#' for membership optimization and \code{"v"} for centroid optimisation.
#'
#' There are seven optimisation algorithms that are currently provided in this
#' package, mainly from the \code{naspaclust}. The optimization algorithm uses
#' the centroid as the parameter to be optimised.
#'
#' Here are the algorithms that can be used:
#' \itemize{
#' \item \code{"classic"} - The classical algorithm of FGWC based on
#' Mason and Jacobson (2007) for centroid optimisation
#' and Runkler and Katz (2006) for membership optimization.
#' \item \code{"abc"} - Optimization using artificial bee colony algorithm
#' based on Karaboga and Basturk (2007) (see also Wijayanto and Purwarianti
#' 2014 and Wijayanto et al. 2016 for FGWC implementation).
#' \item \code{"fpa"} - Optimization using flower pollination algorithm based
#' on (Yang 2012).
#' \item \code{"gsa"} - Optimization using gravitational search algorithm based
#' on Rashedi et al. (2009) and Li and Dong (2017) (see also Pamungkas and
#' Pramana 2019 for FGWC implementation).
#' \item \code{"hho"} - Optimization using harris-hawk optimization with
#' \code{"heidari"} (Heidari et al. 2019) (default), and \code{"bairathi"}
#' (Bairathi and Gopalani 2018).
#' \item \code{"ifa"} - Optimization using intelligent firefly algorithm based
#' on Yang (2009), as well as the intelligent improvement by Fateen and
#' Bonilla-Petriciolet (2013) (see also Nasution et al. 2020 for FGWC
#' implementation).
#' \item \code{"pso"} - Optimization using particle swarm optimization based on
#' Runkler and Katz (2006) and Bansal et al. (2011) for inertia option (see also
#'  Wijayanto and Purwarianti 2014; Putra and Kurniawan 2017; Abdussamad 2020
#'  for FGWC implementation).
#' \item \code{"tlbo"} - Optimization using teaching - learning based
#' optimization based on Rao et al. (2012) and elitism improvement by Rao and
#' Patel (2012).
#' }
#'
#' @importFrom naspaclust fgwcuv
#' @importFrom naspaclust abcfgwc
#'
#' @seealso \code{\link[naspaclust]{fgwcuv}}, \code{\link[naspaclust]{abcfgwc}},
#' \code{\link[naspaclust]{fpafgwc}}, \code{\link[naspaclust]{gsafgwc}},
#' \code{\link[naspaclust]{hhofgwc}}, \code{\link[naspaclust]{ifafgwc}},
#' \code{\link[naspaclust]{psofgwc}}, \code{\link[naspaclust]{tlbofgwc}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname fgwcSTE
#' @references
#' \itemize{
#' \item Abdussamad S (2020). "Evaluation of Implementation Context Based
#' Clustering In Fuzzy Geographically Weighted Clustering-Particle Swarm
#' Optimization Algorithm." Jurnal EECCIS, 14(1), 10--15. ISSN 2460-8122,
#' <https://jurnaleeccis.ub.ac.id/index.php/eeccis/article/view/609>.
#' \item Bairathi D, Gopalani D (2018). "A Novel Swarm Intelligence Based
#' Optimization Method: Harris' Hawk Optimization." In Advances in Intelligent
#' Systems and Computing, 832--842. Springer International Publishing.
#' doi: 10.1007/978-3-030-16660-1_81,
#' <https://doi.org/10.1007/978-3-030-16660-1_81>.
#' \item Bansal JC, Singh PK, Saraswat M, Verma A, Jadon SS, Abraham A (2011).
#' "Inertia Weight strategies in Particle Swarm Optimization." In 2011 Third
#' World Congress on Nature and Biologically Inspired Computing.
#' doi: 10.1109/nabic.2011.6089659,
#' <https://doi.org/10.1109/nabic.2011.6089659>.
#' \item Fateen SK, Bonilla-Petriciolet A (2013). "Intelligent Firefly Algorithm
#'  for Global Optimization." Cuckoo Search and Firefly Algorithm: Theory and
#' Applications, 516, 315--330.
#' \item Heidari AA, Mirjalili S, Faris H, Aljarah I, Mafarja M, Chen H (2019).
#' "Harris hawks optimization: Algorithm and applications." Future Generation
#' Computer Systems, 97, 849--872. doi: 10.1016/j.future.2019.02.028,
#' <https://doi.org/10.1016/j.future.2019.02.028>.
#' \item Karaboga D, Basturk B (2007). "A powerful and efficient algorithm for
#' numerical function optimization: artificial bee colony (ABC) algorithm."
#' Journal of Global Optimization, 39(3), 459--471.
#' doi: 10.1007/s10898-007-9149-x, <https://doi.org/10.1007/s10898-007-9149-x>.
#' \item Li J, Dong N (2017). "Gravitational Search Algorithm with a New
#' Technique." In 2017 13th International Conference on Computational
#' Intelligence and Security (CIS), 516--519. doi: 10.1109/CIS.2017.00120,
#' <https://doi.org/10.1109/CIS.2017.00120>.
#' \item Mason GA, Jacobson RD (2007). "Fuzzy Geographically Weighted
#' Clustering." In Proceedings of the 9th International Conference on
#' Geocomputation, 1--7.
#' \item Nasution BI, Kurniawan R, Siagian TH, Fudholi A (2020). "Revisiting
#' social vulnerability analysis in Indonesia: An optimised spatial fuzzy
#' clustering approach." International Journal of Disaster Risk Reduction, 51,
#' 101801. doi: 10.1016/j.ijdrr.2020.101801,
#' <https://doi.org/10.1016/j.ijdrr.2020.101801>.
#' \item Pamungkas IH, Pramana S (2019). "Improvement Method of Fuzzy
#' Geographically Weighted Clustering using Gravitational Search Algorithm."
#' Journal of Computer Science and Information, 11(1).
#' \item Putra FH, Kurniawan R (2017). "Clustering for Disaster Areas Endemic
#' Dengue Hemorrhagic Fever Based on Factors had Caused in East Java Using Fuzzy
#'  Geographically Weighted Clustering - Particle Swarm Optimization."
#'  Jurnal Aplikasi Statistika & Komputasi Statistik, 8(01), 27. ISSN 2615-1367.
#' \item Rao RV, Patel V (2012). "An elitist teaching-learning-based
#' optimization algorithm for solving complex constrained optimization
#' problems." International Journal of Industrial Engineering Computations,
#' 3(4), 535--560. ISSN 19232926, doi: 10.5267/j.ijiec.2012.03.007,
#' <https://doi.org/10.5267/j.ijiec.2012.03.007>.
#' \item Rao RV, Savsani VJ, Balic J (2012). "Teaching- learning-based
#' optimization algorithm for unconstrained and constrained real-parameter
#' optimization problems." Engineering Optimization, 44(12), 1447--1462.
#' doi: 10.1080/0305215x.2011.652103,
#' <https://doi.org/10.1080/0305215x.2011.652103>.
#' \item Rashedi E, Nezamabadi-pour H, Saryazdi S (2009). "GSA: A Gravitational
#' Search Algorithm." Information Sciences, 179(13).
#' \item Runkler TA, Katz C (2006). "Fuzzy Clustering by Particle Swarm
#' Optimization." In 2006 IEEE International Conference on Fuzzy Systems.
#' doi: 10.1109/fuzzy.2006.1681773,<https://doi.org/10.1109/fuzzy.2006.1681773>.
#' \item Wijayanto AW, Purwarianti A (2014). "Improvement design of fuzzy
#' geo-demographic clustering using Artificial Bee Colony optimization." I
#' n 2014 International Conference on Cyber and IT Service Management (CITSM),
#' 69--74. ISBN 978-1-4799-7975-2.
#' \item Wijayanto AW, Purwarianti A (2014). "Improvement of fuzzy
#' geographically weighted clustering using particle swarm optimization."
#' In 2014 International Conference on Information Technology Systems and
#' Innovation (ICITSI), 7--12. ISBN 978-1-4799-6527-4.
#' \item Wijayanto AW, Purwarianti A, Son LH (2016). "Fuzzy geographically
#' weighted clustering using artificial bee colony: An efficient geo-demographic
#'  analysis algorithm and applications to the analysis of crime behavior in
#'  population." Applied Intelligence, 44(2), 377--398. ISSN 0924-669X.
#' \item Yang X (2014). Nature-Inspired Optimization Algorithms, Elsevier
#' insights. Elsevier Science. ISBN 9780124167452.
#' \item Yang X (2012). "Flower Pollination Algorithm for Global Optimization."
#' In Unconventional Computation and Natural Computation, 240--249. Springer
#' Berlin Heidelberg. doi: 10.1007/978-3-642-32894-7_27,
#' <https://doi.org/10.1007/978-3-642-32894-7_27>.
#' \item Yang X (2009). "Firefly Algorithms for Multimodal Optimization." In
#' Stochastic Algorithms: Foundations and Applications, 169--178. Springer
#' Berlin Heidelberg. doi: 10.1007/978-3-642-04944-6_14,
#' <https://doi.org/10.1007/978-3-642-04944-6_14>.
#' }
#'
#' @export
fgwcSTE <- function(m_sfe,
                    sample_id,
                    data,
                    pop = NULL,
                    distMat = NULL,
                    dMetric = NULL,
                    algorithm = "classic",
                    parameters = NULL) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get distance Matrix
  if (is.null(distMat)) {
    distMat <- getDistMat(sfe, dMetric = dMetric, sample_id = sample_id)
  }

  ## Get population
  if (is.null(pop)) {
    pop <- pop <- rep(1, nrow(data))
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
                      distmat = distMat)

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

  return(fgwc)
}


# ---------------------------------------------------------------------------- #
#  ################# INTERNAL FUNCTIONS ASSOCIATED WITH FGWC ################
# ---------------------------------------------------------------------------- #
#' Internal: Add annotation and geom hex info to fgwc output
#'
#' Adding annotation and/or geom hex info allows easier plotting later.
#'
#' @param fgwc the fgwc output
#' @param sfe the sfe object used as input for `fgwcSTE`
#'
#' @returns description an updated `fgwc$finaldata` data frame.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_addToFGWCout
#'
.int_addToFGWCout <- function(fgwc, sfe, geom_type) {
  annot <- colData(sfe)["annotation"] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Barcode")

  if (geom_type == "hex") {
    geom <- "spotHex"
  }

  geoms <- colGeometry(sfe, geom) %>%
    rownames_to_column(var = "Barcode")

  finaldata <- fgwc$finaldata %>%
    rownames_to_column(var = "Barcode") %>%
    left_join(annot, by = "Barcode") %>%
    left_join(geoms, by = "Barcode") %>%
    column_to_rownames(var = "Barcode")

  return(finaldata)
}


#' Internal Function: Apply Pseudocount
#'
#' This function adds a pseudocount to a numeric vector to prevent zero values.
#'
#' @param values A numeric vector to which the pseudocount will be added.
#' @param pseudocount A numeric value representing the pseudocount to be added
#' (default is \code{1e-12}).
#'
#' @returns A numeric vector with the pseudocount added.
#'
#' @details
#' The pseudocount is used to avoid zero values that can lead to errors in
#' subsequent calculations, such as division by zero or logarithm of zero.
#'
#' @seealso \code{\link[stats]{log}} for logarithmic operations.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal pseudocount
#' @rdname dot-int_applyPseudocount
.int_applyPseudocount <- function(values, pseudocount) {
  values + pseudocount
}


#' Internal Function: Validate Input
#'
#' This function checks if two numeric vectors have the same length.
#'
#' @param x A numeric vector.
#' @param y Another numeric vector.
#'
#' @return Nothing. Stops execution if the input vectors are not of the same
#' length.
#'
#' @details
#' The function ensures that both input vectors are of the same length. If not,
#' it raises an error.
#'
#' @seealso \code{\link[base]{stopifnot}} for input validation.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal validation
#' @rdname dot-int_validateInput
.int_validateInput <- function(x, y) {
  stopifnot(length(x) == length(y))
}


#' Internal Function: Calculate KL Divergence
#'
#' This function calculates the Kullback-Leibler (KL) divergence between two
#' numeric vectors.
#'
#' @param x A numeric vector representing the first probability distribution.
#' @param y A numeric vector representing the second probability distribution.
#' @param pseudocount A numeric value representing the pseudocount to be added
#' to both distributions (default is \code{1e-12}).
#'
#' @return A numeric value representing the KL divergence between the two
#' distributions.
#'
#' @details
#' The function calculates the KL divergence between two probability
#' distributions \code{x} and \code{y}. A pseudocount is added to avoid
#' division by zero or logarithm of zero.
#'
#' @seealso \code{\link[stats]{log}} for logarithmic operations.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal kl-divergence
#' @rdname dot-int_calculateKLdivergence
.int_calculateKLdivergence <- function(x, y, pseudocount = 1e-12) {
  .int_validateInput(x, y)
  x <- .int_applyPseudocount(x, pseudocount)
  y <- .int_applyPseudocount(y, pseudocount)
  sum(x * log(x / y) - x + y)
}


#' Internal Function: Calculate Mean Squared Error (MSE)
#'
#' This function calculates the Mean Squared Error (MSE) between two matrices.
#'
#' @param A.hat A numeric matrix representing the estimated values.
#' @param A A numeric matrix representing the true values.
#'
#' @return A numeric value representing the MSE between the two matrices.
#'
#' @details
#' The function computes the MSE between the estimated matrix (\code{A.hat})
#' and the true matrix (\code{A}).
#'
#' @seealso \code{\link[stats]{mean}} for mean calculations.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal mse
#' @rdname dot-int_calculateMSE
.int_calculateMSE <- function(A.hat, A) {
  if (!inherits(A.hat, "matrix")) {
    A.hat <- as.matrix(A.hat)
  }
  mean((A - A.hat)^2)
}


#' Internal Function: Calculate Mean Kullback-Leibler Divergence (MKL)
#'
#' This function calculates the Mean Kullback-Leibler (MKL) divergence between
#' two matrices.
#'
#' @param A.hat A numeric matrix representing the estimated values.
#' @param A A numeric matrix representing the true values.
#'
#' @return A numeric value representing the MKL divergence between the two matrices.
#'
#' @details
#' The function computes the MKL divergence between the estimated matrix
#' (\code{A.hat}) and the true matrix (\code{A}) by taking the mean KL
#' divergence over all rows.
#'
#' @seealso \code{\link[stats]{mean}} for mean calculations.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal mkl kl-divergence
#' @rdname dot-int_calculateMKL
.int_calculateMKL <- function(A.hat, A) {
  mean(.int_calculateKLdivergence(A.hat, A))
}


#' Internal Function: Initialize Parameters for NMF Optimization
#'
#' This function initializes parameters for Non-negative Matrix Factorization (NMF) optimization.
#'
#' @param A A numeric matrix representing the data to be factorized.
#' @param k_range A numeric vector specifying the range of factor numbers to
#' evaluate.
#' @param loss A character string, either \code{"mse"} (mean squared error) or
#' \code{"mkl"} (mean Kullback-Leibler divergence).
#' @param seed An integer specifying the random seed for reproducibility
#' (default is \code{NULL}).
#'
#' @return Nothing. Initializes the parameters and sets the random seed if
#' provided.
#'
#' @details
#' The function checks that the loss function is valid (\code{"mse"} or
#' \code{"mkl"}), warns about large datasets, sets the random seed if
#' specified, and ensures that the \code{k_range} sequence contains at least
#' three values.
#'
#' @seealso \code{\link[base]{set.seed}} for setting the random seed.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal initialise nmf
#' @rdname dot-int_initialise
.int_initialise <- function(A, k_range, loss, seed) {
  if (!loss %in% c("mse", "mkl")) {
    stop("Invalid loss function. It can only be one of `mse` or `mkl`.")
  }

  if (ncol(A) > 20000) {
    warning("This function can be slow for very large datasets.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(k_range) < 3) {
    stop("k_range sequence must have at least 3 values.")
  }
}


#' Internal Function: Compute Errors for NMF Factor Selection
#'
#' This function computes the reconstruction errors for a given number of
#' factors (\code{k}) using both the input matrix (\code{A}) and a randomized
#' version (\code{A.rand}).
#'
#' @param k An integer specifying the number of factors to use for NMF.
#' @param .A A numeric matrix representing the data to be factorized.
#' @param A.rand A numeric matrix representing the randomized version of the
#' data (\code{A}).
#' @param n_cores An integer indicating the number of cores to use for parallel
#' processing.
#' @param max.iter An integer specifying the maximum number of iterations for
#' NMF.
#' @param loss Character, either \code{"mse"} (mean squared error) or
#' \code{"mkl"} (mean Kullback-Leibler divergence) (default is \code{"mse"}).
#' @param seed An integer specifying the random seed for reproducibility
#'
#' @return A numeric vector with two elements:
#' \describe{
#'   \item{\code{err}}{The reconstruction error for the input matrix
#'   (\code{A}).}
#'   \item{\code{err.rand}}{The reconstruction error for the randomized matrix
#'   (\code{A.rand}).}
#' }
#'
#' @importFrom RcppML nmf
#' @importFrom Matrix Diagonal
#' @details
#' The function uses the \code{nnmf} function from the \code{NNLM} package to
#' perform Non-negative Matrix Factorization (NMF) and calculates the
#' reconstruction errors using either Mean Squared Error (MSE) or Mean
#' Kullback-Leibler (MKL) divergence.
#'
#' @seealso \code{\link[NNLM]{nnmf}} for NMF calculations.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal nmf error-calculation
#' @rdname dot-int_computeErrors
.int_computeErrors <- function(k, A, A.rand, n_cores, max.iter, loss, seed) {
  z <- RcppML::nmf(A = A,
                  k = k,
                  verbose = FALSE,
                  #n.threads = n_cores,
                  #verbose = 0,
                  maxit = max.iter,
                  seed = seed)
  z.rand <- RcppML::nmf(A = A.rand,
                        k = k,
                        verbose = FALSE,
                        #n.threads = n_cores,
                        #verbose = 0,
                        maxit = max.iter,
                        seed = seed)

  A.hat <- with(z, w %*% Matrix::Diagonal(x = d) %*% h)
  A.hat.rand <- with(z.rand, w %*% Matrix::Diagonal(x = d) %*% h)

  if (loss == "mse") {
    err <- .int_calculateMSE(A.hat, A)
    err.rand <- .int_calculateMSE(A.hat.rand, A.rand)
  } else if (loss == "mkl") {
    err <- .int_calculateMKL(A.hat, A)
    err.rand <- .int_calculateMKL(A.hat.rand, A.rand)
  }

  return(c(err, err.rand))
}


#' Internal Function: Find Optimal Number of NMF Factors
#'
#' This function finds the optimal number of factors for Non-negative Matrix
#' Factorization (NMF) by evaluating the reconstruction errors of the input
#' matrix (\code{A}) and a randomized version (\code{A.rand}).
#'
#' @param A A numeric matrix representing the data to be factorized.
#' @param A.rand A numeric matrix representing the randomized version of the
#' data (\code{A}).
#' @param k_range A numeric vector specifying the range of factor numbers to
#' evaluate.
#' @param n_cores An integer indicating the number of cores to use for parallel
#' processing.
#' @param max.iter An integer specifying the maximum number of iterations for
#' NMF.
#' @param loss Character, either \code{"mse"} (mean squared error) or
#' \code{"mkl"} (mean Kullback-Leibler divergence) (default is \code{"mse"}).
#' @param seed An integer specifying the random seed for reproducibility
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{err}}{A numeric vector representing the reduction in
#'   reconstruction error above noise for each factor.}
#'   \item{\code{k}}{An integer representing the optimal number of factors.}
#' }
#'
#' @details
#' The function computes the reconstruction errors for a range of factors
#' (\code{k_range}) using the input matrix and a randomized version, and
#' identifies the optimal number of factors based on the difference in error
#' reduction.
#'
#' @seealso \code{\link[NNLM]{nnmf}} for NMF calculations and
#' \code{\link{.int_computeErrors}} for error computation.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords internal nmf factor-selection
#' @rdname dot-int_findOptimalFactors
.int_findOptimalFactors = function(A,
                                   A.rand,
                                   k_range,
                                   n_cores,
                                   max.iter,
                                   loss,
                                   seed) {
  k.err <- sapply(X = k_range,
                  FUN = .int_computeErrors,
                  A = A,
                  A.rand = A.rand,
                  n_cores = n_cores,
                  max.iter = max.iter,
                  loss = loss,
                  seed = seed)
  # assign("k.err", k.err, envir = .GlobalEnv) # leave thi here for now
  rownames(k.err) <- c("err", "err.rand")
  colnames(k.err) <- k_range

  err.del <- sapply(1:(ncol(k.err) - 1), function(i) k.err[,i] - k.err[,i + 1])
  colnames(err.del) <- colnames(k.err)[2:length(colnames(k.err))]
  err.del.diff <- err.del[1,] - err.del[2,]

  if (sum(err.del.diff < 0) > 0) {
    min.idx <- min(which(err.del.diff < 0))
  } else {
    min.idx <- which.min(err.del.diff) + 1
  }

  return(list(err = err.del.diff, k = k_range[[min.idx]]))
}


#' Internal: Calculate the Fuzzy Partition Coefficient (FPC)
#'
#' Computes the Fuzzy Partition Coefficient (FPC) relative to a fuzzy
#' c-partitioned matrix `u` (essentially the FGWC membership matrix). Measures
#' 'fuzziness' in partitioned clustering.
#'
#' @param u The FGWC membership matrix. A 2D matrix (N, C), where C is the
#' number of clusters and N is the number of data points (locations). Each
#' element represents the membership degree of a data point to a particular
#' cluster.
#'
#' @return A number representing the Fuzzy Partition Coefficient.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords partition coefficient
#' @rdname dot-int_fpcoeff
#'
.int_fpcoeff <- function(u) {
  ## Ensure `u` is a matrix
  u <- as.matrix(u)

  ## Number of locations
  n <- nrow(u)

  ## Calculate fuzzy partition coefficient
  fpc <- sum(diag(t(u) %*% u)) / as.numeric(n)

  return(fpc)
}


#' Internal: Determine the best number of clusters using multiple indexes
#'
#' @param index_data A data frame containing the results of multiple indexes.
#' Columns should include `k`, `FPC`, `CE`, `SC`, `XB`.
#'
#' @return The index_data dataframe with the CompositeScore column added.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords composite score
#' @rdname dot-int_compositeScore
#'
.int_compositeScore <- function(index_data) {
  ## Normalise the indexes (higher is better for FPC and SC, lower is better for CE and XB)
  index_data$FPC <- (index_data$FPC - min(index_data$FPC)) /
    (max(index_data$FPC) - min(index_data$FPC))

  index_data$CE <- (max(index_data$CE) - index_data$CE) /
    (max(index_data$CE) - min(index_data$CE))

  index_data$SC <- (index_data$SC - min(index_data$SC)) /
    (max(index_data$SC) - min(index_data$SC))

  index_data$XB <- (max(index_data$XB) - index_data$XB) /
    (max(index_data$XB) - min(index_data$XB))

  ## Calculate a composite score for each k
  index_data$CompositeScore <- rowMeans(index_data[,c("FPC", "CE", "SC", "XB")])

  ## Delete after
  # assign("index_data", index_data, envir = .GlobalEnv)

  ## Format data frame
  index_data <- data.frame(k = index_data$k,
                           IndexValue = index_data$CompositeScore,
                           IndexType = "Composite")

  return(index_data)
}


#' Internal Function: Find Elbow Point Using Angles
#'
#' This function identifies the elbow point using the angle method.
#'
#' @param index_data A data frame containing columns `k` (number of clusters)
#' and `IndexValue` (metric value).
#'
#' @return An integer representing the optimal number of clusters (`k`).
#'
#' @details
#' This function takes the `index_data` data frame, created internally by a
#' function like \code{fgwc_findOptimumK}, sorts it by `k` (just to ensure it's
#' in the correct order), and then calculates the angle between each pair of
#' adjacent points. The point with the largest angle is considered the elbow
#' point. The function returns the `k` value corresponding to this elbow point.
#'
#' @seealso \code{\link{fgwc_findOptimumK}} for finding the optimal number of
#' clusters.
#'
#' @keywords internal elbow-point
#' @rdname dot-int_findElbowPointAngle
#'
.int_findElbowPointAngle <- function(index_data) {
  ## Ensure the dataframe is sorted by k (just in case)
  index_data <- index_data[order(index_data$k), ]

  ## Calculate the angles between each point and find the maximum angle
  n <- nrow(index_data)
  angles <- rep(NA, n - 2) # n-2 because we have n-2 pairs of k

  for (i in 1:(n - 2)) {
    ## Calculate the angle using the dot product formula for vectors
    v1 <- c(index_data$IndexValue[i + 1] - index_data$IndexValue[i],
            index_data$k[i + 1] - index_data$k[i])
    v2 <- c(index_data$IndexValue[i + 2] - index_data$IndexValue[i + 1],
            index_data$k[i + 2] - index_data$k[i + 1])

    dot_prod <- v1 %*% v2
    mag_v1 <- sqrt(sum(v1 * v1))
    mag_v2 <- sqrt(sum(v2 * v2))

    ## Avoid division by zero
    if (mag_v1 * mag_v2 == 0) {
      angles[i] <- 0
    } else {
      angles[i] <- acos(dot_prod / (mag_v1 * mag_v2))
    }
  }

  ## Find the index of the maximum angle
  elbow_index <- which.max(angles)

  ## The elbow point is the next k after the maximum angle
  elbow_point <- index_data$k[elbow_index + 1]

  return(elbow_point)
}


#' Internal Function: Identify Elbow Point Using the Knee Method
#'
#' This function identifies the elbow point of a curve using the Knee method.
#'
#' @param index_data A data frame containing columns `k` (number of clusters)
#' and `IndexValue` (metric value).
#'
#' @return An integer representing the optimal number of clusters (`k`).
#'
#' @details
#' This function identifies the elbow point by calculating the distance of each
#' point from the line connecting the first and last points. The point with the
#' maximum distance is considered the elbow point.
#'
#' @seealso \code{\link{fgwc_findOptimumK}} for finding the optimal number of
#' clusters.
#'
#' @keywords internal elbow-point
#' @rdname dot-int_findElbowPointKnee
#'
.int_findElbowPointKnee <- function(index_data) {
  ## Ensure the dataframe is sorted by k (just in case)
  index_data <- index_data[order(index_data$k), ]

  ## Export the values
  k_values <- index_data$k
  IndexValue <- index_data$IndexValue

  ## Normalise k and index values
  k_normalised <- (k_values - min(k_values)) / (max(k_values) - min(k_values))
  index_normalised <- (IndexValue - min(IndexValue)) / (max(IndexValue) - min(IndexValue))

  ## Calculate the distances to the line connecting the first and last points
  m <- (index_normalised[length(index_normalised)] - index_normalised[1]) /
    (k_normalised[length(k_normalised)] - k_normalised[1])
  b <- index_normalised[1] - m * k_normalised[1]

  distances <- abs(m * k_normalised - index_normalised + b) / sqrt(m^2 + 1)

  ## Find the index with the maximum distance
  elbow_point_index <- which.max(distances)
  elbow_point <- k_values[elbow_point_index]

  return(elbow_point)
}


# ---------------------------------------------------------------------------- #
#  ##################### NMF FUNCTION COPIED FROM SCATER ####################
# ---------------------------------------------------------------------------- #
## As of 7 Jun 2024, RcppML produces an error in C++. This error is found on
## the latest RcppML update (current GitHub release) but it is not present in
## the previous RcppML version (0.3.7). The scater package must be using the
## latest version. I have installed the RcppML 0.3.7. My use of the package in
## the fgwc_nmfFactorNumber function is not producing the error. I don't know
## what is going on. Until this error is resolved, I will copy and modify scater's
## calculateNMF function.

#' Perform NMF on cell-level data
#'
#' Perform non-negative matrix factorization (NMF) for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param x For \code{calculateNMF}, a numeric matrix of log-expression values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#'
#' For \code{runNMF}, a \linkS4class{SingleCellExperiment} object.
#' @param ... For the \code{calculateNMF} generic, additional arguments to pass to specific methods.
#' For the ANY method, additional arguments to pass to \code{\link[Rtsne]{Rtsne}}.
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#' @return
#' For \code{calculateNMF}, a numeric matrix is returned containing the NMF coordinates for each cell (row) and dimension (column).
#'
#' For \code{runNMF}, a modified \code{x} is returned that contains the NMF coordinates in \code{\link{reducedDim}(x, name)}.
#'
#' In both cases, the matrix will have the attribute \code{"basis"} containing the gene-by-factor basis matrix.
#'
#' @details
#' The function \code{\link[RcppML]{nmf}} is used internally to compute the NMF.
#' Note that the algorithm is not deterministic, so different runs of the function will produce differing results.
#' Users are advised to test multiple random seeds, and then use \code{\link{set.seed}} to set a random seed for replicable results.
#'
#' @seealso
#' \code{\link[RcppML]{nmf}}, for the underlying calculations.
#'
#' \code{\link{plotNMF}}, to quickly visualize the results.
#'
#' @author Aaron Lun
#'

#' @importFrom BiocNeighbors KmknnParam findKNN
#' @importFrom BiocParallel SerialParam
.int_calculate_nmf <- function(x, ncomponents = 2, ntop = 500,
                               subset_row = NULL, scale=FALSE, transposed=FALSE, seed=1, ...)
{
  if (!transposed) {
    x <- scater:::.get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale)
  }
  x <- t(as.matrix(x))

  args <- list(k=ncomponents, verbose=FALSE, seed=seed, ...)
  nmf_out <- do.call(RcppML::nmf, c(list(x), args))

  # RcppML doesn't use transposed data
  nmf_x <- t(nmf_out$h)
  rownames(nmf_x) <- colnames(x)
  colnames(nmf_x) <- paste0("NMF", seq_len(ncol(nmf_x)))
  nmf_basis <- nmf_out$w
  rownames(nmf_basis) <- rownames(x)
  colnames(nmf_basis) <- paste0("NMF", seq_len(ncol(nmf_basis)))
  attr(nmf_x, "basis") <- nmf_basis

  nmf_x
}




