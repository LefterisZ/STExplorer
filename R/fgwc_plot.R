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
#' @rdname plotFGWC_singleMap
#' @aliases plotFGWC_single
#'
#' @export
plotFGWC_singleMap <- function(fgwc,
                               m_sfe,
                               sample_id,
                               colours = NULL) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Prepare data to plot
  data <- .int_fgwcPlotDataMap(fgwc = fgwc, sfe = sfe, mode = "single")

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
            aes(geometry = geometry,
                fill = as.factor(Cluster)),
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
#' @rdname plotFGWC_multiMap
#' @aliases plotFGWC_multi
#'
#' @export
plotFGWC_multiMap <- function(fgwc,
                              m_sfe,
                              sample_id,
                              palette = "YlGnBu") {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Prepare data to plot
  data <- .int_fgwcPlotDataMap(fgwc = fgwc, sfe = sfe, mode = "multi") %>%
    tidyr::pivot_longer(cols = !.data$geometry,
                        names_to = "cluster",
                        values_to = "membership")

  ggplot() +
    geom_sf(data = data,
            aes(geometry = geometry,
                fill = membership),
            colour = "grey30",
            show.legend = TRUE) +
    scale_fill_distiller(palette = palette, limits = c(0, 1)) +
    facet_wrap(~cluster) +
    labs(fill = "Cluster\nmembership") +
    theme_void() +
    theme(legend.position = "right")

}


#' Plot Membership percentages from FGWC per annotation or winning cluster
#'
#' This function plots the results of Fuzzy Geographically Weighted Clustering
#' (FGWC) for multiple clusters. It plots the membership percentage of each
#' cluster in each location as violin plots.
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
#' @rdname plotFGWC_multiViolin
#' @aliases plotFGWC_multi
#'
#' @export
plotFGWC_multiViolin <- function(fgwc) {

  ## Prepare data to plot
  data <- .int_fgwcPlotDataViolin(fgwc = fgwc, sfe = sfe) %>%
    tidyr::pivot_longer(cols = -all_of(c("annotation", "Cluster")),
                        names_to = "cluster",
                        values_to = "membership")

  ggplot(data = data,
         aes(x = cluster, y = membership)) +
    geom_violin(aes(fill = cluster),
                width = 1,
                show.legend = TRUE,
                trim = FALSE) +
    geom_boxplot(width = 0.1, colour = "black", alpha = 0.2,
                 outlier.size = 0.5, outlier.alpha = 0.2) +
    scale_fill_viridis_d() +
    facet_wrap(~annotation) +
    labs(x = "",
         y = "Membership %",
         fill = "Cluster") +
    theme_classic() +
    theme(legend.position = "right",
          strip.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 17))

}


#' Plot NMF Factors as a Map
#'
#' This function plots the Non-Negative Matrix Factorization (NMF) factors as a map.
#'
#' @param nmf The NMF result object. As generated by STExplorer's `fgwc_nmf`
#' function.
#' @param m_sfe The SpatialFeatureExperiment (SFE) or metaSFE object.
#' @param sample_id The sample ID.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' plotFGWC_nmfFactorsMap(nmf_result, m_sfe, sample_id)
#' }
#'
#' @seealso \code{\link{nmf}}, \code{\link{plotFGWC_single}},
#' \code{\link{plotFGWC_multi}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @keywords nmf factorization spatial ggplot2 visualization
#' @family plotting functions
#' @rdname plotFGWC_nmfFactorsMap
#' @aliases plotFGWC_nmfFactorsMap
#'
#' @export
plotFGWC_nmfFactorsMap <- function(nmf, m_sfe, sample_id) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Prepare data to plot
  data <- .int_nmfPlotDataMap(nmf = nmf, sfe = sfe)

  ## Plot
  ggplot(data = data) +
    geom_sf(aes(geometry = geometry, fill = score)) +
    facet_wrap(~ Factors) +
    scale_fill_viridis_c() +
    labs(fill = "Factor\nScore") +
    theme_void()
}


plotFGWC_nmfFactorsPie <- function() {}


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
  rownames(annot_row) <- make.unique(markers$gene.name) # add gene names in rownames
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
    geom_sf(aes(geometry = geometry, fill = gTruth),
            alpha = 0.8) +
    geom_sf(aes(geometry = geometry_subC, colour = subclust),
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
#' @rdname dot-int_fgwcPlotDataMap
#'
.int_fgwcPlotDataMap <- function(fgwc, sfe, mode = c("single", "multi")) {
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
  } else if (mode == "both") {

  } else {
    stop("\n'mode' argument accepts only 'single' or 'multi' as values.\n")
  }

  return(fgwc_membership)
}


#' Internal Function: Generate Data for FGWC Violin Plot
#'
#' This function generates data suitable for plotting FFGWC results as violins.
#'
#' @param fgwc An object containing FGWC results, typically obtained from the
#' `fgwcSTE` function.
#' @param sfe A SpatialFeatureExperiment object containing spatial coordinates
#' and features.
#'
#' @return A data frame containing information for generating FGWC plots. In
#' single-cluster mode, the data frame includes cluster information for each
#' spot. In multi-cluster mode, the entire data frame is returned.
#'
#' @details The function generates data suitable for plotting FGWC results as
#' violin plots.
#'
#' @seealso \code{\link[SpatialFeatureExperiment]{colGeometry}},
#' \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords FGWC, clustering, spatial coordinates
#'
#' @rdname dot-int_fgwcPlotDataViolin
#'
.int_fgwcPlotDataViolin <- function(fgwc, sfe) {
  clusts.number <- 1:ncol(fgwc$membership)
  fgwc_membership <- data.frame(fgwc$membership)
  colnames(fgwc_membership) <- paste0("Cluster_", clusts.number)
  rownames(fgwc_membership) <- colnames(sfe)
  fgwc_membership <- data.frame(fgwc_membership,
                                Cluster = fgwc$cluster)
  if (!is.null(colData(sfe)[["annotation"]])) {
    fgwc_membership$annotation = colData(sfe)$annotation
  } else {
    fgwc_membership <- fgwc_membership %>%
      mutate(annotation = Cluster)
  }

  return(fgwc_membership)
}


#' Internal Function: generate data for NMF factor scores
#'
#' This function generates data suitable for plotting FFGWC results as violins.
#'
#' @param nmf An object containing NMF results, typically obtained from the
#' `fgwc_nmf` function.
#' @param sfe A SpatialFeatureExperiment object containing spatial coordinates
#' and features.
#'
#' @return A data frame containing information for generating NMF plots.
#'
#' @details The function generates data suitable for plotting NMF results as
#' map.
#'
#' @seealso \code{\link[SpatialFeatureExperiment]{colGeometry}},
#' \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords FGWC, clustering, spatial coordinates
#'
#' @rdname dot-int_nmfPlotDataMap
#'
.int_nmfPlotDataMap <- function(nmf, sfe) {
  scores <- (nmf) %>%
    as.data.frame() %>%
    rownames_to_column()

  geoms <- colGeometry(sfe, "spotHex") %>%
    rownames_to_column()

  # Combine spatial coordinates with scores
  scores <- dplyr::left_join(geoms, scores, by = "rowname") %>%
    column_to_rownames() %>%
    pivot_longer(cols = -geometry, names_to = "Factors", values_to = "score")

  return(scores)
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
                     by = dplyr::join_by("ensg.ID")) %>%
    ## Remove ENSG.IDs
    dplyr::select(-"ensg.ID") %>%
    ## Bring gene names column to the front
    dplyr::relocate("gene.name") %>%
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
