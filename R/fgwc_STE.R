#' Identify the optimal number of factors
#'
#' This function uses different statistical methods to select the optimal
#' number of factors for a given dataset.
#'
#' @param name description
#'
#' @return An integer. The number of optimal factors.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords clustering spatial-expression nmf
#' @family clustering functions
#' @rdname fgwc_nmfFactorNumber
#' @aliases fgwc_nmfFactorNumber
#'
#' @export
fgwc_nmfFactorNumber <- function() {}


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
#' @importFrom scater calculateNMF
#' @importFrom dplyr select
#'
#' @export
fgwc_nmf <- function(m_sfe,
                     sample_id,
                     assay = "logcounts",
                     top_hvgs,
                     ncomponents = 2,
                     ntop = 600,
                     subset_row = NULL,
                     scale = FALSE,
                     ...) {
  ## SFE or metaSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get data
  nmf_input <- assay(sfe, assay)[rownames(sfe) %in% top_hvgs,]

  ## Run NMF
  nmf <- scater::calculateNMF(x = nmf_input,
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

#' #' Fuzzy Geographicaly Weighted Clustering
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
#' a data frame with 2 columns.
#' @param pop an n*1 vector contains population.
#' @param dMetric Character string specifying the distance metric.
#' @param distMat an n*n distance matrix between regions.
#' @param algorithm algorithm used for FGWC
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
#'
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


  return(fgwc)
}
