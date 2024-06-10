# Prepare lists ----
samples <- names(msfe@sfe_data)

best_k_nmf <- list()
sfe_nmf_list <- list()
best_k_fgwc <- list()
fgwc_param_list <- list()
fgwc_list <- list()

marker_heatmap_list <- list()
subHeatmap_list <- list()

# Find optimum factor number ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Find optimum number of Factors
  result <- fgwc_nmfFactorNumber(m_sfe = msfe,
                                 sample_id = s,
                                 assay = "logcounts",
                                 top_hvgs = top_hvgs[[s]],
                                 k_range = seq(2, 10, 1),
                                 n_cores = 1,
                                 do_plot = FALSE,
                                 seed = 1,
                                 loss = "mse",
                                 max.iter = 250)

  best_k_nmf[[s]] <- result

  print(plotFGWC_factorSelection(result))

  ## Housekeeping
  rm(result)
}

# Calculate NMF ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Run NMF
  sfe_nmf_list[[s]] <- fgwc_nmf(m_sfe = msfe,
                                sample_id = s,
                                top_hvgs = top_hvgs[[s]],
                                ncomponents = best_k_nmf[[s]][["k"]])
}

# Find best k for clustering ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Find best number of clusters
  fgwc_param_list[[s]] <- fgwc_params(algorithm = "classic", ncluster = 5)

  best_k_fgwc[[s]] <- fgwc_findOptimumK(fgwc_in = sfe_nmf_list[[s]],
                                        k_range = 2:10,
                                        index_type = "FPC",
                                        elbow_method = "knee",
                                        m_sfe = msfe,
                                        sample_id = s,
                                        algorithm = "classic",
                                        parameters = fgwc_param_list[[s]])

  ## update the parameters input
  fgwc_param_list[[s]] <- fgwc_params(algorithm = "classic",
                                      ncluster = best_k_fgwc[[s]])
}

# Run FGWC ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Run FGWC
  fgwc_list[[s]] <- fgwcSTE(m_sfe = msfe,
                            sample_id = s,
                            data = sfe_nmf_list[[s]],
                            dMetric = "euclidean",
                            algorithm = "classic",
                            parameters = fgwc_param_list[[s]])
}

# Plot single clusters ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Plot single clusters
  p1 <- plotFGWC_singleMap(fgwc = fgwc_list[[s]],
                           m_sfe = msfe,
                           sample_id = s)

  p2 <- plotQC_spotsAnnotation(msfe,
                               sample_id = s,
                               type = "hex")

  print(p1 + p2)

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_single_NMF-", ncomp, "_clust-",
                         ncluster, "_a-", a, "_m-", m, ".svg"),
         device = "svg",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)
}

# Plot multi clusters ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Plot multi clusters
  p1 <- plotFGWC_multiMap(fgwc = fgwc_list[[s]],
                          m_sfe = msfe,
                          sample_id = s)

  print(p1)

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_multi_NMF-", ncomp, "_clust-",
                         ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot memberships in violins ----
for (s in samples) {
  ## Plot memberships in violins
  p1 <- plotFGWC_multiViolin(fgwc = fgwc_list[[s]],
                             m_sfe = msfe,
                             sample_id = s)

  print(p1)

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_multiViolin_NMF-", ncomp, "_clust-",
                         ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot FGWC pie-doughnuts ----
for (s in samples) {
  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m

  ## Matching FGWC clusters and annotation of locations
  plotFGWC_pie(fgwc = fgwc_list[[s]],
               m_sfe = msfe,
               sample_id = s,
               mapping = aes(pie = cluster, donut = annotation))
  ggplot2::ggsave(paste0(folder, s, "_pieClustAnn_NMF-", ncomp, "_clust-",
                         ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)

  ## Matching FGWC clusters and NMF factors
  plotFGWC_pie(fgwc = fgwc_list[[s]],
               m_sfe = msfe,
               sample_id = s,
               mapping = aes(pie = cluster, donut = factors))
  ggplot2::ggsave(paste0(folder, s, "_pieClustFact_NMF-", ncomp, "_clust-",
                         ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot map of factor scores ----
for (s in samples) {
  ## Plot a map of the NMF factor scores for each location
  p1 <- plotFGWC_nmfFactorsMap(nmf = sfe_nmf_list[[s]],
                               m_sfe = msfe,
                               sample_id = s)

  print(p1)

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_factorsMap_NMF-", ncomp, "_clust-",
                         ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot heatmap of factor scores ----
for (s in samples) {
  ## Plot a factor heatmap alongside spot annotations and/or FGWC clusters
  plotFGWC_nmfFactorsHeatmap(fgwc = fgwc_list[[s]],
                             loc_annot = "both",
                             order_rows = "cluster")

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_factorsHeatMap_NMF-", ncomp,
                         "_clust-", ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot heatmap of Metagene scores ----
for (s in samples) {
  ## Plot a map of the NMF factor scores for each location
  plotFGWC_nmfMetagenesHeatmap(fgwc = fgwc_list[[s]])

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_metaGenesHeatMap_NMF-", ncomp,
                         "_clust-", ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot marker gene heatmap ----
for (s in samples) {
  ## Load the cell markers as a named list
  # markers <- list(sample1 = markers_sample1,
  #                 sample2 = markers_sample2)

  cluster_no <- 2

  ## Plot the heatmap
  heatmap <- plotFGWC_markersHeatmap(fgwc = fgwc_list[[s]],
                                     m_sfe = msfe,
                                     sample_id = s,
                                     markers = markers[[s]],
                                     cluster_no = cluster_no,
                                     cutree_cols = 5)

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_markersHeatMap", cluster_no,
                         "_NMF-", ncomp, "_clust-", ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot subcluster map ----
for (s in samples) {
  clust <- 3
  ## Plot the subclusters
  plotFGWC_subClust(heatmap = marker_heatmap_list[[s]],
                    k = 5,
                    clust = clust,
                    m_sfe = msfe,
                    sample_id = s)

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_subClust", clust,
                         "_NMF-", ncomp, "_clust-", ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# Plot subcluster heatmap ----
for (s in samples) {
  clust <- 3
  cluster_no <- 4
  ## Only plot the sub-cluster heatmap
  plotFGWC_subHeatmap(heatmap = marker_heatmap_list[[s]],
                      k = 5,
                      markers = markers[[s]],
                      m_sfe = msfe,
                      sample_id = s,
                      cluster_no = cluster_no)

  ## Save the heatmap as an object
  subHeatmap_list[[s]] <- plotFGWC_subHeatmap(heatmap = marker_heatmap_list[[s]],
                                              k = 5,
                                              markers = markers[[s]],
                                              m_sfe = msfe,
                                              sample_id = s,
                                              cluster_no = cluster_no)

  folder <- paste0("./data/graphics_out/benchmarking/FGWC_in-STExplorer/clustering/", s, "/")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
  ggplot2::ggsave(paste0(folder, s, "_subClustHeatmap_Clust", clust, "subClust", cluster_no,
                         "_NMF-", ncomp, "_clust-", ncluster, "_a-", a, "_m-", m, ".svg"),
                  device = "svg",
                  width = grDevices::dev.size(units = "in")[1],
                  height = grDevices::dev.size(units = "in")[2],
                  units = "in",
                  dpi = 300)
}

# ---------------------------------------------------------------------------- #
#' Plot FGWC NMF Factors Heatmap
#'
#' This function generates a clustered heatmap of factor scores obtained from
#' FGWC (Fuzzy Geographically Weighted Clustering) analysis. It allows for
#' optional annotations based on clusters or other metadata.
#'
#' @param fgwc An object of class `fgwc` containing the results of FGWC analysis.
#' @param loc_annot A character string specifying the type of annotations to include
#'   in the heatmap. Possible values are:
#'   \itemize{
#'     \item `"annotation"`: Include only the annotation column.
#'     \item `"cluster"`: Include only the cluster column.
#'     \item `"both"`: Include both annotation and cluster columns.
#'     \item `"none"`: Include no annotations.
#'   }
#'   Defaults to `"both"`.
#'
#' @return A heatmap plot of factor scores with optional annotations.
#' @importFrom dplyr select contains all_of any_of mutate
#' @importFrom pheatmap pheatmap
#' @importFrom cols4all c4a
#'
#' @examples
#' \dontrun{
#' # Assuming `fgwc_object` is a valid FGWC object
#' plotFGWC_nmfFactorsHeatmap(fgwc_object, loc_annot = "both")
#' plotFGWC_nmfFactorsHeatmap(fgwc_object, loc_annot = "annotation")
#' plotFGWC_nmfFactorsHeatmap(fgwc_object, loc_annot = "none")
#' }
#'
#' @export
plotFGWC_nmfFactorsHeatmap <- function(fgwc,
                                       loc_annot = c("annotation", "cluster",
                                                     "both", "none")) {
  ## Check fgwc argument is of class fgwc
  .int_checkFGWCClass(fgwc)

  ## Prepare data for heatmap
  data_in <- .int_getFactorData(fgwc = fgwc)

  ## Prepare location annotations
  if (loc_annot != "none") {
    annotations_and_colours <- .int_getAnnotsAndColours(fgwc = fgwc,
                                                        loc_annot = loc_annot)
  } else {
    annotations_and_colours <- list(annotations = NA, colours = NA)
  }

  ## Create heatmap of Factor scores
  pheatmap(data_in,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_row = annotations_and_colours$annotations,
           annotation_colors = annotations_and_colours$colours,
           fontsize_row = 3)
}



#' Internal: Check if an Object is of Class 'fgwc'
#'
#' This function checks if a given object is of class 'fgwc'.
#'
#' @param object The object to check.
#' @return A logical value indicating whether the object is of class 'fgwc'.
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_checkFGWCClass
#'
.int_checkFGWCClass <- function(object) {
  if (!inherits(object, "fgwc")) {
    stop("The `fgwc` argument must an object of class `fgwc` as generated by ",
         "the `fgwc_STE` function.")
  }
}


#' Internal: Extract Factor Data
#'
#' This function extracts factor data from FGWC output
#'
#' @param fgwc fgwc class object as generated by `fgwcSTE`
#'
#' @returns a data frame of NMF factors
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_getFactorData
#'
.int_getFactorData <- function(fgwc) {
  fgwc$finaldata %>%
    dplyr::select(contains("Factor"))
}


#' Internal: Get Annotations and Colours
#'
#' This function returns annotations label and colours for the heatmap
#'
#' @param fgwc gwc class object as generated by `fgwcSTE`
#' @param loc-annot which annotations to add
#'
#' @returns a list
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_getAnnotsAndColours
#'
.int_getAnnotsAndColours <- function(fgwc, loc_annot) {
  annotations <- fgwc$finaldata %>%
    dplyr::select(dplyr::any_of(c("cluster", "annotation"))) %>%
    dplyr::mutate(cluster = as.factor(cluster))

  col_clust <- length(unique(annotations$cluster))
  annot_colours <- list(cluster = getColours(col_clust))
  names(annot_colours$cluster) <- unique(annotations$cluster)

  if ("annotation" %in% colnames(annotations)) {
    col_annot <- length(unique(annotations$annotation))
    annot_colours[["annotation"]] <- c4a(palette = "carto.pastel", n = col_annot)
    names(annot_colours$annotation) <- unique(annotations$annotation)
  }

  if (loc_annot %in% c("cluster", "annotation")) {
    annotations <- annotations %>%
      dplyr::select(all_of(loc_annot))
    annot_colours <- annot_colours[[loc_annot]]
  }

  list(annotations = annotations, colours = annot_colours)
}


#' @importFrom BiocNeighbors KmknnParam findKNN
#' @importFrom BiocParallel SerialParam
.int_calculate_nmf <- function(x, ncomponents = 2, ntop = 500,
                           subset_row = NULL, scale=FALSE, transposed=FALSE, seed=1, ...)
{
  if (!transposed) {
    x <- scater:::.get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale)
  }
  x <- t(as.matrix(x))

  # args <- list(k=ncomponents, verbose=FALSE, seed=seed)
  # nmf_out <- do.call(RcppML::nmf, c(list(x), args))

  nmf_out <- RcppML::nmf(A = x,
                         k = ncomponents,
                         verbose = FALSE,
                         #n.threads = n_cores,
                         #verbose = 0,
                         maxit = 250,
                         seed = seed)

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


# ---------------------------------------------------------------------------- #
nmf_input <- assay(sfe, "logcounts")[rownames(sfe) %in% top_hvgs[["V1_2"]],]
