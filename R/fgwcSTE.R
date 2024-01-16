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
#' @return A named vector of FGWC parameters or optimization algorithm
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
                        ncluster,
                        kind = "u",
                        m = 2,
                        distance = "manhattan",
                        order = 1,
                        alpha = 0.5,
                        a = 1,
                        b = 1,
                        max.iter = 500,
                        error = 1e-5,
                        randomN = 1,
                        uij = NA,
                        vi = NA) {
  ## Check algorithm arguments
  algorithm <- match.arg(algorithm)

  if (algorithm == "classic") {
    ## A named vector that consists of FGWC parameter, or a vector that consists
    ##   of optimization algorithm parameter.
    params <- c(kind = kind,
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
  else if (algorithm != "classic") {
    stop("\nCurrently this function supports only algortihm type 'classic'.\n",
         "For the other optimisation algorithms, refer to the vignette with\n",
         "the default parameters for any optimisation algortihm or refer to\n",
         "the help page (?fgwc_params) to find the optimisation algorithm's\n",
         "help under section 'See also'.")
  }

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
#' @param fgwc_param a vector that consists of FGWC parameter
#' (see \code{\link{fgwcuv}} for parameter details)
#' @param opt_param a vector that consists of optimization algorithm parameter
#' (see \code{\link{fgwcuv}} for parameter details)
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
#' @importFrom naspaclust fgwc
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
                    fgwc_param = NULL,
                    opt_param = NULL) {
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
  if (is.null(fgwc_param)) {
    fgwc_param <- c(kind = 'u', ncluster = 2, m = 2,
                    distance = 'euclidean', order = 2,
                    alpha = 0.7, a = 1, b = 1,
                    max.iter = 500, error = 1e-5, randomN = 1)
  }
  if (is.null(opt_param)) {
    opt_param <- c(vi.dist = 'uniform',
                   npar = 10, par.no = 2, par.dist = 'euclidean', par.order = 2,
                   pso = TRUE, same = 10, type = 'sim.annealing',
                   ei.distr = 'normal', vmax = 0.7, wmax = 0.9, wmin = 0.4,
                   chaos = 4, x0 = 'F', map = 0.7, ind = 1, skew = 0, sca = 1)
  }

  ## Run FGWC
  fgwc <- naspaclust::fgwc(data = data,
                           pop = pop,
                           distmat = distMat,
                           algorithm = algorithm,
                           fgwc_param = fgwc_param,
                           opt_param = opt_param)

  return(fgwc)
}

#' Plot Single Cluster Result for Fuzzy Geographically Weighted Clustering
#' (FGWC)
#'
#' This function plots the result of Fuzzy Geographically Weighted Clustering
#' (FGWC) for a single cluster. The selected cluster is the cluster with the
#' highest membership percentage for each location.
#'
#' @param fgwc The result of FGWC clustering.
#' @param m_sfe The spatial feature experiment (SFE) or metaSFE object.
#' @param sample_id The sample ID.
#' @param colours A vector of colours for the clusters. If NULL, default
#' colours will be used.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' nmf <- fgwc_nmf(m_sfe, sample_id, top_hvgs, ncomponents = 2,
#' ntop = 600)
#' fgwc_result <- fgwc(m_sfe, sample_id, data = nmf)
#' plotFGWC_single(fgwc_result, m_sfe, sample_id)
#' }
#'
#' @seealso
#' \code{\link{fgwc_nmf}}, \code{\link{plotFGWC_multi}},
#' \code{\link{plotFGWC_heatmap}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords clustering fuzzy-logic spatial ggplot2 visualization fgwc
#' @family plotting functions
#' @rdname plotFGWC_single
#' @aliases plotFGWC_single
#'
#' @export
plotFGWC_single <- function(fgwc,
                            m_sfe,
                            sample_id,
                            colours = NULL) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Prepare data to plot
  data <- .int_fgwcPlotData(fgwc = fgwc, sfe = sfe, mode = "single")

  ## Fetch colours
  if (is.null(colours)) {
    col.No = length(unique(fgwc$cluster))
    annot_cols <- getColours(col.No)
  } else {
    annot_cols <- colours
  }

  ## Plot single cluster
  ggplot() +
    geom_sf(data = data,
            aes(geometry = data$geometry,
                fill = as.factor(data$Cluster)),
            colour = "grey30",
            show.legend = TRUE) +
    scale_fill_manual(values = annot_cols) +
    labs(title = NULL,
         fill = "Cluster") +
    theme_void() +
    theme(legend.position = "right")

}

#' Plot Multiple Cluster Results for Fuzzy Geographically Weighted Clustering
#' (FGWC)
#'
#' This function plots the results of Fuzzy Geographically Weighted Clustering
#' (FGWC) for multiple clusters. It plots the membership percentage of each
#' cluster in each location.
#'
#' @param fgwc The result of FGWC clustering.
#' @param m_sfe The spatial feature experiment (SFE) or metaSFE object.
#' @param sample_id The sample ID.
#' @param palette The color palette for the cluster memberships.
#' Default is "YlGnBu".
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' fgwc_result <- fgwc_nmf(m_sfe, sample_id, top_hvgs, ncomponents = 2,
#' ntop = 600)
#' plotFGWC_multi(fgwc_result, m_sfe, sample_id)
#' }
#'
#' @seealso
#' \code{\link{fgwc_nmf}}, \code{\link{plotFGWC_single}},
#' \code{\link{plotFGWC_heatmap}}
#'
#' @details
#' The function produces a panel of maps. One map per cluster. The colours in
#' all maps are scaled to be the same. Meaning that the colour for 50% in map 1
#' is going to be the same with colour for 50% in map 3.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 scale_fill_distiller
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords clustering fuzzy-logic spatial ggplot2 visualization fgwc
#' @family plotting functions
#' @rdname plotFGWC_multi
#' @aliases plotFGWC_multi
#'
#' @export
plotFGWC_multi <- function(fgwc,
                           m_sfe,
                           sample_id,
                           palette = "YlGnBu") {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Prepare data to plot
  data <- .int_fgwcPlotData(fgwc = fgwc, sfe = sfe, mode = "multi") %>%
    tidyr::pivot_longer(cols = !.data$geometry,
                        names_to = "cluster",
                        values_to = "membership")

  ggplot() +
    geom_sf(data = data,
            aes(geometry = .data$geometry,
                fill = .data$membership),
            colour = "grey30",
            show.legend = TRUE) +
    scale_fill_distiller(palette = palette, limits = c(0, 1)) +
    facet_wrap(~cluster) +
    labs(fill = "Cluster\nmembership") +
    theme_void() +
    theme(legend.position = "right")

}


#' Plot Heatmap for Fuzzy Geographically Weighted Clustering (FGWC)
#'
#' This function plots a heatmap for the gene expression of marker genes for a
#' specific cluster.
#'
#' @param fgwc The result of FGWC clustering.
#' @param m_sfe The spatial feature experiment (SFE) or metaSFE object.
#' @param sample_id The sample ID.
#' @param markers A data frame with 4 columns: gene.name, ensg.ID, Type,
#' and Subtype.
#' @param cluster_no The cluster number for which the heatmap is generated.
#' @param cutree_cols Optional, cutree result for columns.
#'
#' @return A pheatmap object.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' fgwc_result <- fgwc_nmf(m_sfe, sample_id, top_hvgs, ncomponents = 2,
#' ntop = 600)
#'
#' markers_data <- read.csv("markers_data.csv") # Provide your markers data
#' cluster_no <- 1 # Provide the cluster number for which heatmap is desired
#'
#' plotFGWC_heatmap(fgwc_result, m_sfe, sample_id, markers_data, cluster_no)
#' }
#'
#' @seealso
#' \code{\link{fgwc_nmf}}, \code{\link{plotFGWC_single}},
#' \code{\link{plotFGWC_multi}}
#'
#' @details
#' Additional details about the function or its behaviour can be added here.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords clustering fuzzy-logic spatial heatmap visualization fgwc
#' @family plotting functions
#' @rdname plotFGWC_heatmap
#' @aliases plotFGWC_heatmap
#'
#' @importFrom grDevices colorRampPalette
#'
#' @export
plotFGWC_heatmap <- function(fgwc,
                             m_sfe,
                             sample_id,
                             markers,
                             cluster_no,
                             cutree_cols = NA) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## 'markers' argument should be a data frame with 4 columns: gene.name,
  ##  ensg.ID, Type, and Subtype
  if (length(markers) == 0) {
    stop("The 'markers' argument is empty.",
         " Please provide a 4-column data frame.",
         "\n If in doubt what the data frame should include, look at the ",
         "format section of the markers example dataset typing ?markers")
  }

  ## Prepare expression data for heatmap
  markers <- .int_markers(markers = markers)
  marker_counts <- .int_fgwcMarkerCounts(sfe = sfe, markers = markers)

  ## Add clusters to gene counts
  marker_clusts <- data.frame("cluster" = fgwc$cluster,
                              marker_counts)

  ## Prepare pheatmap input
  pheat_in <- .int_pheatInput(marker_clusts = marker_clusts,
                              cluster_no = cluster_no,
                              markers = markers)

  ## Heatmap
  ## Annotate rows
  annot_row <- data.frame(Type = markers$Type, # get annotations
                          Subtype = markers$Subtype)
  rownames(annot_row) <- markers$gene.name # add gene names in rownames
  ## Remove genes that might not be present in the input
  annot_row <- annot_row %>%
    filter(rownames(.) %in% rownames(pheat_in)) %>%
    arrange(.data[["Type"]])
  col_type <- length(unique(markers$Type))
  col_subT <- length(unique(markers$Subtype))
  annot_colours <- list(Type = getColours(col_type),
                        Subtype = c4a(palette = "carto.pastel", n = col_subT))
  names(annot_colours$Type) <- unique(markers$Type)
  names(annot_colours$Subtype) <-  unique(markers$Subtype)

  ## group rows
  gaps_row <- utils::head(as.numeric(cumsum(table(annot_row$Type))), -1)

  ## Heatmap colour and set it around zero
  paletteLength <- 25
  c4a_palette <- c4a(palette = "tol.sunset")[c(1,5,9)]
  ph.color <- colorRampPalette(c(c4a_palette[1],
                                 "white",
                                 c4a_palette[3]))(paletteLength)
  ph.breaks <- c(seq(min(pheat_in),
                     0,
                     length.out = ceiling(paletteLength/2) + 1),
                 seq(max(pheat_in)/paletteLength,
                     max(pheat_in),
                     length.out = floor(paletteLength/2)))

  ## Plot heatmap
  heatmap <- pheatmap(pheat_in,
                      color = ph.color,
                      breaks = ph.breaks,
                      annotation_colors = annot_colours,
                      annotation_row = annot_row,
                      gaps_row = gaps_row,
                      show_colnames = FALSE,
                      cluster_rows = FALSE,
                      cluster_cols = TRUE,
                      cutree_cols = cutree_cols,
                      main = paste0("Cluster ", cluster_no))

  return(heatmap)
}


#' Plot FGWC Sub-Clusters in Map
#'
#' This function generates a map of sub-clusters based on the fuzzy
#' geographically weighted clustering (FGWC) results.
#'
#' @param heatmap The FGWC heatmap object.
#' @param k The number of sub-clusters.
#' @param clust The main cluster number.
#' @param m_sfe Either a SpatialfeatureExperiment or a
#' SpatialfeatureExperiment object containing spatial transcriptomics data.
#' @param sample_id The sample.
#' @param annot_cols Vector of colours for ground truth annotations. If NULL,
#' colours are assigned automatically.
#' @param subClust_cols Vector of colours for sub-clusters. If NULL, colours are
#'  generated automatically.
#'
#' @return A ggplot object visualizing the sub-clusters in the spatial context.
#'
#' @details
#' This function takes the results of the FGWC clustering and produces a spatial
#'  map of sub-clusters for a specific FGWC cluster. It overlays the identified
#'  sub-clusters on top of the original spatial features, allowing users to
#'  visualise the distribution of sub-clusters within the main cluster areas.
#'
#' @importFrom sf st_polygon
#' @importFrom ggplot2 scale_colour_manual
#'
#' @examples
#' \dontrun{
#' # Plot sub-clusters with default colors.
#' plotFGWC_subClust(heatmap, k = 5, clust = 3, m_sfe, sample_id = "JBO019")
#'
#' # Plot sub-clusters with custom colors.
#' plotFGWC_subClust(heatmap, k = 5, clust = 3, m_sfe, sample_id = "JBO019",
#'                   annot_cols = c("red", "blue"),
#'                   subClust_cols = c("orange","green","yellow","red","blue"))
#' }
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @export
plotFGWC_subClust <- function(heatmap, k, clust,
                              m_sfe, sample_id,
                              annot_cols = NULL, subClust_cols = NULL) {
  ## Plot sub-clusters in map

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Extracting sub-clusters
  subclusts <- .int_getSubClusts(heatmap = heatmap, k = k)

  ## Extracting necessary information
  sfe_data <- colData(sfe)
  spot_hex_geometry <- colGeometry(sfe, "spotHex")

  ## Creating sub-cluster map
  subclust_map <- data.frame(
    Barcode = rownames(sfe_data),
    geometry = spot_hex_geometry,
    gTruth = sfe_data$annotation,
    geometry_subC = spot_hex_geometry
    ) %>%
    dplyr::left_join(subclusts, by = "Barcode") %>%
    dplyr::mutate(
      subclust = ifelse(.data$subclust > 0,
                        paste0(clust, LETTERS[.data$subclust]),
                        NA)
      ) %>%
    rename("geometry_subC" = "geometry.1")
  subclust_map$geometry_subC[is.na(subclust_map$subclust)] <- st_polygon()

  ## Extract colours
  if (is.null(annot_cols)) {
    colNo1 <- .int_getLabelsNo(subclust_map$gTruth)
    colour.annot <- getColours(colNo1)
  } else {
    colour.annot <- annot_cols
  }
  if (is.null(subClust_cols)) {
    colNo2 <- .int_getLabelsNo(subclust_map$subclust)
    colour.subClust <- .int_getColSubset(colNo1, colNo2)
  } else {
    colour.subClust <- subClust_cols
  }

  ggplot(data = subclust_map) +
      geom_sf(aes(geometry = .data$geometry, fill = .data$gTruth),
              alpha = 0.8) +
      geom_sf(aes(geometry = .data$geometry_subC, colour = .data$subclust),
              fill = NA,
              linewidth = 1.1) +
      scale_fill_manual(values = colour.annot,
                        na.value = "grey95") +
      scale_colour_manual(values = colour.subClust,
                          na.value = "grey95") +
      labs(title = paste0("Cluster ", clust),
           fill = "Ground\nTruth",
           colour = "Sub-clusters") +
      theme_void()
}


#' Plot an FGWC Sub-Cluster Heatmap of gene expression
#'
#' This function generates a heatmap for the gene expression profiles within a
#' specific sub-cluster identified by fuzzy geographically weighted clustering
#' (FGWC).
#'
#' @param heatmap The FGWC heatmap object.
#' @param k The number of sub-clusters used in FGWC.
#' @param markers A data frame with gene markers information, including columns
#' "gene.name", "ensg.ID"  "Type", and "Subtype".
#' @param m_sfe Either a a SpatialfeatureExperiment or a
#' SpatialfeatureExperiment object containing spatial transcriptomics data.
#' @param sample_id The sample ID used for FGWC.
#' @param cluster_no The specific sub-cluster number for which the heatmap is
#' generated.
#' @param cutree_cols Number of columns for cutree to cut the sample tree of the
#'  heatmap. If NA, the tree will not be cut.
#'
#' @return A pheatmap object representing the gene expression heatmap for the
#' specified sub-cluster.
#'
#' @details
#' This function takes the results of the FGWC clustering, extracts the gene
#' expression profiles within the specified sub-cluster, and generates a
#' heatmap. The heatmap provides insights into the expression patterns of
#' selected markers within the identified sub-cluster.
#'
#' @examples
#' \dontrun{
#' # Plot a heatmap for sub-cluster 3 using default colors.
#' plotFGWC_subHeatmap(heatmap, k = 5, markers, m_sfe, sample_id = "JBO019",
#' cluster_no = 3)
#'
#' # Plot a heatmap for sub-cluster 2 with custom cutree colors.
#' plotFGWC_subHeatmap(heatmap, k = 5, markers, m_sfe, sample_id = "JBO019",
#' cluster_no = 2, cutree_cols = 4)
#' }
#'
#' @export
plotFGWC_subHeatmap <- function(heatmap,
                                k,
                                markers,
                                m_sfe,
                                sample_id,
                                cluster_no,
                                cutree_cols = NA) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Extracting sub-clusters
  subclusts <- .int_getSubClusts(heatmap = heatmap, k = k)

  ## Prepare expression data for heatmap
  markers <- .int_markers(markers = markers)
  marker_counts <- .int_fgwcMarkerCounts(sfe = sfe, markers = markers) %>%
    rownames_to_column(var = "Barcode")

  ## Add clusters to gene counts
  marker_clusts <- left_join(subclusts, marker_counts, by = "Barcode") %>%
    column_to_rownames(var = "Barcode")

  ## Prepare pheatmap input
  pheat_in <- .int_pheatInput(marker_clusts = marker_clusts,
                              cluster_no = cluster_no,
                              markers = markers)

  ## Heatmap
  ## Annotate rows
  annot_row <- data.frame(Type = markers$Type, # get annotations
                          Subtype = markers$Subtype)
  rownames(annot_row) <- markers$gene.name # add gene names in rownames
  ## Remove genes that might not be present in the input
  annot_row <- annot_row %>%
    filter(rownames(.) %in% rownames(pheat_in)) %>%
    arrange(.data$Type)
  col_type <- length(unique(markers$Type))
  col_subT <- length(unique(markers$Subtype))
  annot_colours <- list(Type = getColours(col_type),
                        Subtype = c4a(palette = "carto.pastel", n = col_subT))
  names(annot_colours$Type) <- unique(markers$Type)
  names(annot_colours$Subtype) <-  unique(markers$Subtype)

  ## group rows
  gaps_row <- utils::head(as.numeric(cumsum(table(annot_row$Type))), -1)

  ## Heatmap colour and set it around zero
  paletteLength <- 25
  c4a_palette <- c4a(palette = "tol.sunset")[c(1,5,9)]
  ph.color <- colorRampPalette(c(c4a_palette[1],
                                 "white",
                                 c4a_palette[3]))(paletteLength)
  ph.breaks <- c(seq(min(pheat_in),
                     0,
                     length.out = ceiling(paletteLength/2) + 1),
                 seq(max(pheat_in)/paletteLength,
                     max(pheat_in),
                     length.out = floor(paletteLength/2)))

  ## Plot heatmap
  heatmap <- pheatmap(pheat_in,
                      color = ph.color,
                      breaks = ph.breaks,
                      annotation_colors = annot_colours,
                      annotation_row = annot_row,
                      gaps_row = gaps_row,
                      show_colnames = FALSE,
                      cluster_rows = FALSE,
                      cluster_cols = TRUE,
                      cutree_cols = cutree_cols,
                      main = paste0("Sub-Cluster ", cluster_no))

  return(heatmap)
}


# ---------------------------------------------------------------------------- #
#  ############## INTERNAL FUNCTIONS ASSOCIATED WITH THE FGWC ###############
# ---------------------------------------------------------------------------- #
#' Internal Function: Generate Data for FGWC Plot
#'
#' This function generates data suitable for plotting FFGWC results.
#'
#' @param fgwc An object containing FGWC results, typically obtained from the
#' `fgwcSTE` function.
#' @param sfe A SpatialFeatureExperiment object containing spatial coordinates
#' and features.
#' @param mode Character, indicating the mode for generating plot data. Options
#' are "single" (default) for a single-cluster mode or "multi" for a
#' multi-cluster mode.
#'
#' @return A data frame containing information for generating FGWC plots. In
#' single-cluster mode, the data frame includes cluster information for each
#' spot. In multi-cluster mode, the entire data frame is returned.
#'
#' @details The function generates data suitable for plotting FGWC results. In
#' single-cluster mode, the data frame includes cluster information for each
#' spot. In multi-cluster mode, the entire data frame is returned. The mode is
#' specified using the 'mode' argument.
#'
#' @seealso \code{\link[SpatialFeatureExperiment]{colGeometry}},
#' \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords FGWC, clustering, spatial coordinates
#'
#' @rdname dot-int_fgwcPlotData
#'
.int_fgwcPlotData <- function(fgwc, sfe, mode = c("single", "multi")) {
  clusts.number <- 1:ncol(fgwc$membership)
  fgwc_membership <- data.frame(fgwc$membership)
  colnames(fgwc_membership) <- paste0("Cluster_", clusts.number)
  rownames(fgwc_membership) <- colnames(sfe)
  fgwc_membership <- data.frame(fgwc_membership,
                                geometry = colGeometry(sfe, "spotHex"))
  if (mode == "single") {
    fgwc_membership$Cluster <- fgwc$cluster
  } else if (mode == "multi") {
    return(fgwc_membership)
  } else {
    stop("\n'mode' argument accepts only 'single' or 'multi' as values.\n")
  }

  return(fgwc_membership)
}


#' Internal Function: Retrieve Normalised Counts for gene Markers
#'
#' This function retrieves normalised counts for a specified set of gene
#' markers from a SpatialFeatureExperiment object.
#'
#' @param sfe A SpatialFeatureExperiment object containing normalised counts.
#' @param markers A data frame containing information about FGWC markers,
#' including gene names and corresponding ENSG.IDs.
#'
#' @return A data frame with normalised counts for the specified FGWC markers.
#'
#' @details The function extracts normalised counts for a set of FGWC markers
#' from the provided SpatialFeatureExperiment object. It filters the counts
#' based on the ENSG.IDs of the specified markers, rearranges columns for
#' better readability, and provides a data frame with gene names as row names
#' and observations as variables.
#'
#' @seealso \code{\link[SpatialFeatureExperiment]{colData}},
#' \code{\link{as.data.frame}}, \code{\link[dplyr]{filter}},
#' \code{\link[dplyr]{left_join}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords FGWC, normalised counts, markers, SpatialFeatureExperiment
#'
#' @rdname dot-int_fgwcMarkerCounts
#'
.int_fgwcMarkerCounts <- function(sfe, markers) {
  ## Get normalised counts
  counts <- assay(sfe, "logcounts") %>%
    as(., "matrix") %>%
    as.data.frame() %>%
    ## Get Genes Of Interest (GOI)
    dplyr::filter(rownames(.) %in% markers$ensg.ID) %>%
    ## Get GOI ENSG.IDs to a column
    tibble::rownames_to_column(var = "ensg.ID") %>%
    ## Join with annotations
    dplyr::left_join(.,
              markers[,c("gene.name", "ensg.ID")],
              by = dplyr::join_by(.data$ensg.ID)) %>%
    ## Remove ENSG.IDs
    dplyr::select(-.data$ensg.ID) %>%
    ## Bring gene names column to the front
    dplyr::relocate(.data$gene.name) %>%
    ## Populate rownames with gene names
    tibble::column_to_rownames("gene.name") %>%
    ## Transpose to make it observation x variable
    t() %>%
    as.data.frame()

  return(counts)
}


#' Internal Function: Prepare Input Data for Plotting a Heatmap with Feature
#' Grouping
#'
#' This function prepares input data for plotting a heatmap with feature
#' grouping based on cluster or subcluster results.
#'
#' @param marker_clusts A data frame containing information about marker genes
#' and their cluster assignments.
#' @param cluster_no Numeric, the cluster or subcluster number for which to
#' prepare the input data.
#' @param markers A data frame containing information about FGWC markers,
#' including gene names.
#'
#' @return A scaled and ordered data frame suitable for generating a heatmap
#' with cell type.
#'
#' @details The function extracts information about marker genes and their
#' cluster or subcluster assignments from the provided data frame. It scales
#' the data, sets NaN values to zero, orders the data based on FGWC marker
#' genes, and drops any introduced NAs, resulting in a data frame ready for
#' plotting a heatmap with cell type.
#'
#' @seealso \code{\link{scale}}, \code{\link[dplyr]{filter}},
#' \code{\link[dplyr]{select}}, \code{\link[tidyr]{drop_na}}
#'
#' @importFrom tidyr drop_na
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords heatmap, cell type, clustering
#'
#' @rdname dot-int_pheatInput
#'
.int_pheatInput <- function(marker_clusts, cluster_no, markers) {
  ## Check for clusters or subclusters
  if ("cluster" %in% colnames(marker_clusts)) {
    input <- marker_clusts %>%
      dplyr::filter(.data$cluster == cluster_no) %>%
      dplyr::select(-all_of("cluster"))
  } else if ("subclust" %in% colnames(marker_clusts)) {
    input <- marker_clusts %>%
      dplyr::filter(.data$subclust == cluster_no) %>%
      dplyr::select(-all_of("subclust"))
  }

  ## Proceed with scaling and further wrangling
  input <- input %>%
    scale() %>%
    t() %>%
    as.data.frame()

  ## Set NaN values into zeros
  input[is.na(input)] <- 0

  ## Order and drop NAs introduced
  input <- input %>%
    .[markers$gene.name,] %>%
    tidyr::drop_na()

  return(input)
}


#' Process cell type Markers
#'
#' This function processes cell type markers, ensuring required column
#' names, arranging by marker types, and removing duplicates.
#'
#' @param markers A data frame containing information about cell type
#' markers, including gene names and types.
#'
#' @return A processed data frame with validated column names, arranged by
#' marker types, and duplicates removed.
#'
#' @details The function validates the column names of the input data frame,
#' ensuring the presence of "gene.name," "ensg.ID," "Type," and "Subtype"
#' columns. It arranges the markers by types and removes duplicate entries
#' based on ENSG.IDs, resulting in a processed data frame ready for further
#' analysis.
#'
#' @seealso \code{\link{.int_validateColumnNames}},
#' \code{\link{.int_arrangeByType}}, \code{\link{.int_removeDuplicates}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords cell type, markers, validation, duplicates
#'
#' @rdname dot-int_markers
#'
.int_markers <- function(markers) {
  .int_validateColumnNames(markers,
                           c("gene.name", "ensg.ID", "Type", "Subtype"))

  markers <- .int_arrangeByType(markers)
  markers <- .int_removeDuplicates(markers, "ensg.ID")

  return(markers)
}


#' Validate Column Names in a Data Frame
#'
#' This function validates the column names in a data frame, ensuring they
#' match the expected column names.
#'
#' @param df A data frame to validate.
#' @param expected_columns A character vector specifying the expected column
#' names.
#'
#' @return NULL if validation is successful; otherwise, an error is raised.
#'
#' @details The function checks if the actual column names in the provided data
#' frame match the expected column names. If there is a mismatch, it raises an
#' error indicating the expected and actual column names for user correction.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords validation, data frame, column names
#'
#' @rdname dot-int_validateColumnNames
#'
.int_validateColumnNames <- function(df, expected_columns) {
  actual_columns <- colnames(df)
  if (!all(actual_columns %in% expected_columns) ||
      length(actual_columns) != length(expected_columns)) {
    stop("Invalid column names. Expected: ",
         paste(expected_columns, collapse = ", "),
         "; Actual: ",
         paste(actual_columns, collapse = ", "))
  }
}


#' Arrange cell type Markers by Type
#'
#' This function arranges cell type markers by their types in
#' ascending order.
#'
#' @param markers A data frame containing information about cell type
#' markers, including gene names and types.
#'
#' @return A data frame with markers arranged by their types in ascending order.
#'
#' @details The function uses the `dplyr::arrange` function to sort the markers
#' based on the cell types in ascending order, resulting in a data frame where
#' markers are organised by their types for improved readability.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords cell type, markers, arrangement, types
#'
#' @rdname dot-int_arrangeByType
#'
.int_arrangeByType <- function(markers) {
  markers <- markers %>%
    dplyr::arrange(.data$Type)

  return(markers)
}


#' Remove Duplicates from cell type Markers
#'
#' This function removes duplicate entries from cell type markers based
#' on a specified column.
#'
#' @param markers A data frame containing information about cell type
#' markers, including gene names and types.
#' @param column The column based on which duplicates are identified
#' and removed.
#'
#' @return A data frame with duplicate entries removed based on the specified
#' column.
#'
#' @details The function uses the `BiocGenerics::duplicated` function to
#' identify and tag duplicate entries in the specified column. It then filters
#' out those entries, resulting in a data frame with unique cell type
#' markers based on the chosen column.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords cell type, markers, duplicates, removal
#'
#' @importFrom BiocGenerics duplicated
#'
#' @rdname dot-int_removeDuplicates
#'
.int_removeDuplicates <- function(markers, column) {
  markers <- markers %>%
    dplyr::mutate(dups = BiocGenerics::duplicated({{ column }})) %>%
    dplyr::filter(!.data$dups)

  return(markers)
}


#' Get the Number of Unique Labels in a Vector
#'
#' This function calculates the number of unique labels in a vector, excluding
#' missing values.
#'
#' @param vector A vector for which to calculate the number of unique labels.
#'
#' @return An integer representing the number of unique labels in the vector.
#'
#' @details The function utilises the `unique` function to extract unique l
#' abels from the vector. It then counts the non-missing labels and returns the
#' total number.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords unique labels, vector, counting
#'
#' @rdname dot-int_getLabelsNo
#'
.int_getLabelsNo <- function(vector) {
  labels <- unique(vector)
  out <- sum(!is.na(labels))

  return(out)
}


#' Get a Subset of Colours from a Palette
#'
#' This function retrieves a subset of colours from a palette based on
#' specified indices.
#'
#' @param no1 The starting index for the subset.
#' @param no2 The ending index for the subset.
#'
#' @return A vector of colours representing the subset from the palette.
#'
#' @details The function uses the `getColours` function to obtain a palette of
#' colours. It then extracts a subset of colours starting from the specified
#' index `no2` up to the end of the palette.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords subset, colours, palette
#'
#' @rdname dot-int_getColSubset
#'
.int_getColSubset <- function(no1, no2) {
  colours <- getColours(no1 + no2)
  out <- colours[no2:length(colours)]

  return(out)
}

#' Get Subclusters from Heatmap
#'
#' Assigns subclusters to barcodes based on the input heatmap and the specified
#' number of clusters.
#'
#' @param heatmap Data containing clustering information.
#' @param k Number of subclusters.
#'
#' @return A data frame with barcode and corresponding subcluster assignment.
#'
.int_getSubClusts <- function(heatmap, k) {
  # Get the subclusters
  subclusts <- stats::cutree(heatmap$tree_col, k = k) %>%
    .[heatmap$tree_col$order] %>%
    as.data.frame() %>%
    dplyr::rename("subclust" = ".") %>%
    tibble::rownames_to_column(var = "Barcode")

  return(subclusts)
}
