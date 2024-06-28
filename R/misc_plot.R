#' Visualise Neighbour Graphs and Associated Images
#'
#' This function generates a grid of neighbour graphs for spatial
#' transcriptomics data, optionally overlaying image data.
#'
#' @param msfe The spatial transcriptomics data as a
#' MetaSpatialFeatureExperiment object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param res Character vector specifying the resolution of the images to be
#' used ('lowres', 'hires', 'fullres').
#' @param plotImage Logical. Whether to overlay image data on the neighbour
#' graphs.
#' @param ... Additional parameters passed to internal functions.
#'
#' @details The function creates a grid of neighbour graphs for spatial
#' transcriptomics data. Each graph represents the spatial relationships
#' between spots. The optional overlay of image data provides context to the
#' spatial distribution.
#'
#' @return A grid of neighbour graphs, with or without overlaid image data.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial transcriptomics, neighbour graph, grid, image overlay
#'
#' @rdname plotNeighbourGraph
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(msfe_example)
#' plotNeighbourGraph(msfe_example, sample_id = "JBO019",
#'                    res = 'lowres', plotImage = TRUE)
#' }
#' @export
plotNeighbourGraph <- function(msfe,
                               sample_id = TRUE,
                               res = c("lowres", "hires", "fullres"),
                               plotImage = TRUE,
                               ...) {
  ## Check arguments
  res <- match.arg(res)

  ## Select samples
  ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)

  ## Prepare neighbour graph data
  nb_data_list <-  lapply(msfe@sfe_data[ids], .int_getNBdata)
  ## Fetch image data and transform to raster
  if (plotImage) {
    image_list <- lapply(msfe@sfe_data[ids], .int_getImgDtMSFE, image_id = res)
  }
  ## Get capture area limits
  limits_list <- lapply(msfe@sfe_data[ids], .int_getImgLims)

  ## Plot
  plots <- lapply(ids,
                  .int_plotNBgraph,
                  plotImage = plotImage,
                  image_list = image_list,
                  nb_data_list = nb_data_list,
                  limits_list = limits_list)

  ## Determine the number of columns in the grid (can be customized)
  num_cols <- ceiling(sqrt(length(plots)))

  ## Combine the individual plots into a grid
  gridExtra::grid.arrange(grobs = plots, ncol = num_cols)


}


#' Map Gene Expression
#'
#' A function to plot gene expression data on spatial feature
#' experiment (SFE) objects as a map.
#'
#' @param m_sfe An object of class 'SpatialFeatureExperiment' or
#'              'MetaSpatialFeatureExperiment' containing spatial
#'              data.
#' @param genes A vector of gene names or ENSG_IDs to be plotted.
#' @param sample_id A character string specifying the sample ID for
#'                  analysis. If NULL the first sample is selected.
#' @param assay Character vector specifying the type of assay data
#'              to plot. Options are "counts" for raw counts and
#'              "logcounts" for log2-normalized counts. Default is
#'              "counts".
#' @param minmax Numeric vector of length 2 specifying the minimum
#'               and maximum expression values to include in the
#'               plot. Default is c(0, Inf), meaning all values
#'               above 0 will be included.
#' @param type Character vector specifying the type of spatial data
#'             to use for plotting. Options are "spot" for spot
#'             geometry and "hex" for hexagon geometry. Default is
#'             "spot".
#' @param res Character vector specifying the resolution of the image
#'            to fetch for plotting. Options are "lowres" for
#'            low-resolution, "hires" for high-resolution,
#'            "fullres" for full-resolution, and "none" to skip
#'            fetching image data. Default is "none".
#' @param fill_args Additional arguments for customizing the fill
#'                  colour palette.
#' @param alpha Numeric value specifying the transparency level of
#'              the plotted polygons. Default is 0.3.
#' @param ... Additional arguments to be passed to the plotting
#'            functions.
#'
#' @importFrom purrr reduce
#' @importFrom tidyterra geom_spatraster_rgb
#' @importFrom ggplot2 geom_sf labs coord_sf facet_wrap theme_void
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' plotGeneExpression(m_sfe = sfe_object,
#'                    genes = c("Gene1", "Gene2"),
#'                    sample_id = "sample123",
#'                    assay = "logcounts",
#'                    minmax = c(0, 100),
#'                    type = "hex",
#'                    res = "hires",
#'                    fill_args = list(option = "magma",
#'                                     na.value = "grey"),
#'                    alpha = 0.5)
#' }
#'
#' @seealso \code{\link{addSpatialData}}, \code{\link{plotSpatialData}}
#'
#' @return A plot of gene expression data on spatial feature
#'         experiment objects.
#'
#' @author Eleftherios Zormpas
#'
#' @export
plotGeneExpression <- function(m_sfe,
                               genes,
                               sample_id = NULL,
                               assay = c("counts", "logcounts"),
                               minmax = c(0, Inf),
                               type = c("spot", "hex", "cntd"),
                               res = c("lowres", "hires", "fullres", "none"),
                               fill_args = list(),
                               alpha = 0.3,
                               ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
    type <- "spotHex"
  } else if (type == "spot") {
    type <- "spotPoly"
  } else if (type == "cntd") {
    type <- "spotCntd"
  }

  ## Fetch image if needed
  if (!res == "none") {
    ## Fetch image data and transform to raster
    image <- .int_getImgDt(sfe = sfe, sample_id = sample_id, image_id = res)
    ## Get capture area limits
    limits_list <- .int_getImgLims(sfe = sfe)
  }

  ## Create legend title
  if (assay == "counts") {
    fill <- "Raw counts"
  } else if (assay == "logcounts") {
    fill <- "Log2-Normalised\ncounts"
  }

  ## Set fill arguments if not provided
  if (isEmpty(fill_args)) {
    fill_args <- list(option = "viridis",
                      na.value = "grey")
  }

  ## Fetch data
 data <- .int_prepareDataMap(sfe = sfe,
                             assay = assay,
                             genes = genes,
                             minmax = minmax,
                             type = type)

  ## Plot
 p <- ggplot()

 if (!res == "none") {
   p <- p +
     tidyterra::geom_spatraster_rgb(data = image[[1]])
 }

 if (type == "cntd") {
   p <- p +
     ggplot2::geom_sf(data = data,
                      aes(geometry = geometry,
                          fill = expression),
                      alpha = alpha) +
     do.call(scale_fill_viridis_c, c(list(), fill_args)) +
     ggplot2::labs(fill = fill)
 } else {
   p <- p +
     ggplot2::geom_sf(data = data,
                      aes(geometry = geometry,
                          colour = expression),
                      alpha = alpha,
                      size = 0.5) +
     do.call(scale_colour_viridis_c, c(list(), fill_args)) +
     ggplot2::labs(colour = fill)
 }

 p <- p +
   ggplot2::coord_sf() +
   ggplot2::facet_wrap(~gene, scales = "fixed", ncol = 1) +
   ggplot2::theme_void()

 p
}


#' Plot a Heatmap of Gene Expression Data
#'
#' This function plots a heatmap using the `pheatmap` package, providing
#' options to subset and annotate the data.
#'
#' @param m_sfe The \code{SpatialFeatureExperiment} or
#' \code{MetaSpatialFeatureExperiment} object containing spatial expression
#' data.
#' @param sample_id A character string indicating the name of the sample. It is
#' required when the `m_sfe` argument is provided with a MetaSFE object.
#' Defaults to `NULL` which selects the first sample from the `m_sfe` object.
#' @param subset_col subset cols (locations) using a boolean vector of length
#' equal to the number of locations.
#' @param subset_row subset rows (genes) using a boolean vector of length
#' equal to the number of genes.
#' @param loc_annot a named list of annotations for the columns (locations). If
#' left empty, the function looks if there is an `annotation` column
#' in the `colData`. If there is not, then, no annotations are added.
#' @param order_col order the columns (locations) based on the alphabetical
#' order of the values provided in the `loc_annot` argument. The value needs to
#' be a character string of the name of the annotation to be used for ordering.
#' This name MUST match the names in the named list provided in the `loc_annot`
#' argument. Defaults to `NULL` which leads to columns being clustered
#' by `pheatmap`.
#' @param scale character indicating if the values should be centred and scaled
#' in either the row direction or the column direction, or none. Corresponding
#' values are "row", "column" and "none". (this is a `pheatmap` parameter)
#' @param ... Additional arguments passed to `pheatmap`.
#'
#' @returns Plots a heatmap using `pheatmap`. It can be saved as a `pheatmap`
#' object.
#'
#' @details
#' The function assumes that the annotations provided in the `loc_annot`
#' argument are in the same order as the locations within the `SFE` object.
#' Please make sure this is the case.
#'
#' The named list provided in the `loc_annot` argument MUST have the values as
#' vectors (i.e., `list(annot = vector_of_annotations)`). The
#' `vector_of_annotations` in the above example MUST have the same length as
#' the number of locations in the `SFE` object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname plotHeatmap
#'
#' @examples
#' \dontrun{
#' # Assuming `sfe` is a valid SFE object
#' # Assuming there is also cluster information from running `fgwc_STE`
#'
#' # Preparing the annotations
#' loc_annots <- list(cluster = fgwc$cluster,
#'                    annotation = colData(sfe)$annotation)
#' # The values provided in the list MUST be vectors.
#'
#' # Preparing to order the columns based on cluster number
#' order_cols <- "cluster"
#'
#' plotHeatmap(sfe,
#'             loc_annot = loc_annots,
#'             order_col = order_cols)
#' }
#'
#' @export
plotHeatmap <- function(m_sfe,
                        sample_id = NULL,
                        assay = "logcounts",
                        subset_col = NULL,
                        subset_row = NULL,
                        loc_annot = list(),
                        order_col = NULL,
                        scale = "row",
                        ...) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Prepare heatmap data
  data_in <- .int_getGeneHeatmapData(sfe = sfe,
                                     assay = assay,
                                     subset_row = subset_row,
                                     subset_col = subset_col)

  ## Prepare location annotations
  if (loc_annot != "none") {
    annots_and_colours <- .int_getAnnotsAndColoursSFE(sfe = sfe,
                                                      loc_annot = loc_annot,
                                                      order_col = order_col)
  } else {
    annotations_and_colours <- list(annotations = NA, colours = NA)
  }

  ## Order data input by annotation or cluster
  if (is.null(order_col)) {
    cluster_cols = TRUE
  } else {
    data_in <- data_in[,rownames(annotations_and_colours$annotations)]
    cluster_cols = FALSE
  }

  ## Add labels
  show_colnames <- ifelse(ncol(data_in) <= 20, TRUE, FALSE)
  show_rownames <- ifelse(nrow(data_in) <= 50, TRUE, FALSE)

  ## Add custom colours
  colour <- rev(cols4all::c4a(palette = "brewer.br_bg", n = 9))

  ## Create heatmap of Factor scores
  pheatmap(data_in,
           color = colour,
           scale = scale,
           cluster_rows = TRUE,
           cluster_cols = cluster_cols,
           show_rownames = show_rownames,
           show_colnames = show_colnames,
           annotation_col = annots_and_colours$annotations,
           annotation_colors = annots_and_colours$colours,
           fontsize_row = 3,
           ...)
}
# ---------------------------------------------------------------------------- #
#  ############# INTERNAL FUNCTIONS ASSOCIATED WITH MISC PLOTS ##############
# ---------------------------------------------------------------------------- #
#' Internal Function: Fetch Image Data and Transform to Raster
#'
#' This internal function retrieves image data for a given
#' SpatialFeatureExperiment and transforms it into a raster format.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param image_id Character vector specifying the image ID.
#'
#' @return A raster object containing the image data.
#'
#' @details
#' This function is an internal utility to fetch image data from a
#' SpatialFeatureExperiment using the specified image ID. The retrieved image
#' is then transformed into a raster format for further processing.
#'
#' @seealso
#' \code{\link[SpatialFeatureExperiment]{getImg}},
#' \code{\link[SpatialFeatureExperiment]{SFE-image}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial image raster
#'
#' @rdname dot-int_getImgDtMSFE
#'
#' @aliases .int_getImgDtMSFE
#'
#' @importFrom SpatialFeatureExperiment imgRaster
#' @importFrom SpatialFeatureExperiment getImg
#'
.int_getImgDtMSFE <- function(sfe, image_id) {
  ## Fetch image data and transform to raster
  image <- getImg(sfe, image_id = image_id) %>% imgRaster()

  return(image)
}


#' Internal Function: Extract Neighbour Graph Data
#'
#' This internal function extracts neighbour graph data from a
#' SpatialFeatureExperiment.
#'
#' @param sfe A SpatialFeatureExperiment object.
#'
#' @return A Simple Features (sf) object representing the neighbour graph data.
#'
#' @details
#' This function is an internal utility to extract neighbour graph data from a
#' SpatialFeatureExperiment. It utilizes the colGraph function to obtain the
#' neighbours and spatial coordinates, and then converts the data to a Simple
#' Features (sf) object.
#'
#' @seealso
#' \code{\link{colGraph}}, \code{\link{spatialCoords}}, \code{\link{nb2lines}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial neighbour graph sf
#'
#' @rdname dot-int_getNBdata
#' @aliases .int_getNBdata
#' @importFrom spdep nb2lines
#'
.int_getNBdata <- function(sfe) {
  nb_data <-  as(nb2lines(colGraph(sfe)$neighbours,
                          coords = spatialCoords(sfe)),
                 "sf")

  return(nb_data)
}


#' Internal Function: Plot Neighbour Graph with Optional Tissue Image
#'
#' This internal function plots the neighbour graph with an optional overlay
#' of tissue image for a selected sample.
#'
#' @param id The sample ID.
#' @param plotImage Logical, indicating whether to plot the tissue image.
#' Default is TRUE.
#' @param image_list List containing image data for each sample.
#' @param nb_data_list List containing neighbour graph data for each sample.
#' @param limits_list List containing capture area limits for each sample.
#'
#' @details
#' This function is an internal utility to visualize the neighbour graph for a
#' selected sample. It allows the option to overlay the tissue image on the
#' graph. The capture area limits are considered for proper spatial alignment.
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot}}, \code{\link[tidyterra]{geom_spatraster_rgb}},
#' \code{\link[ggplot2]{geom_sf}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial visualization plotting ggplot2
#'
#' @rdname dot-int_plotNBgraph
#' @aliases .int_plotNBgraph
#' @importFrom tidyterra geom_spatraster_rgb
#'
.int_plotNBgraph <- function(id,
                             plotImage,
                             image_list,
                             nb_data_list,
                             limits_list) {
  ## Plot with or without the tissue image
  if (plotImage) {
    p <- ggplot() +
      tidyterra::geom_spatraster_rgb(data = image_list[[id]]) +
      ggplot2::geom_sf(data = nb_data_list[[id]])
  } else {
    p <- ggplot() +
      ggplot2::geom_sf(data = nb_data_list[[id]])
  }

  p + ggplot2::labs(subtitle = id) +
    ggplot2::lims(x = limits_list[[id]][[1]],
                  y = limits_list[[id]][[2]]) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    theme(plot.subtitle = ggplot2::element_text(hjust = 0.5))

  # return(p)
}


#' Internal: prepare Data Map
#'
#' A function to prepare data for plotting gene expression on spatial
#' feature experiment (SFE) objects.
#'
#' @param sfe An object of class 'SpatialFeatureExperiment' containing spatial
#'            data.
#' @param assay Character vector specifying the type of assay data
#'              to use. Options are "counts" for raw counts and
#'              "logcounts" for log2-normalized counts.
#' @param genes A vector of gene names or ENSG_IDs to be plotted.
#' @param minmax Numeric vector of length 2 specifying the minimum
#'               and maximum expression values to include in the
#'               plot.
#' @param type Character vector specifying the type of spatial data
#'             to use. Options are "spot" for spot geometry and "hex"
#'             for hexagon geometry.
#'
#' @return A data frame containing gene expression data prepared for
#'         plotting on spatial feature experiment objects.
#'
#' @importFrom purrr reduce
#' @importFrom dplyr select filter
#' @importFrom tidyr pivot_longer
#' @importFrom Matrix t
#'
.int_prepareDataMap <- function(sfe, assay, genes, minmax, type) {
  # genes_df <- assay(sfe, assay) %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   dplyr::select(tidyr::all_of(genes)) %>%
  #   tibble::rownames_to_column(var = "rownames")
  genes_df <- assay(sfe, assay) %>%
    Matrix::t() %>%
    as.data.frame()

  genes_lst <- lapply(genes, FUN = function(x) {
    genes_df %>%
      dplyr::select(tidyr::all_of(x)) %>%
      tibble::rownames_to_column(var = "rownames") %>%
      filter(.data[[x]] > minmax[1] & .data[[x]] < minmax[2])
  })

  geoms <- data.frame(geometry = colGeometry(sfe, type)) %>%
    tibble::rownames_to_column(var = "rownames")

  geoms_lst <- list()
  geoms_lst[[1]] <- geoms

  data_lst <- c(geoms_lst, genes_lst)

  out <- purrr::reduce(data_lst, dplyr::left_join, by = 'rownames') %>%
    tibble::column_to_rownames(var = "rownames") %>%
    tidyr::pivot_longer(cols = -"geometry",
                        names_to = "gene",
                        values_to = "expression")

  return(out)
}


#' Internal: Extract Membership Data
#'
#' This function extracts membership data from FGWC output
#'
#' @param fgwc fgwc class object as generated by `fgwcSTE`
#' @param subset_col subset cols (locations) using a boolean vector of length
#' equal to the number of locations.
#' @param subset_row subset rows (genes) using a boolean vector of length
#' equal to the number of genes.
#'
#' @returns a data frame of cluster membership %
#'
#' @importFrom tidyr drop_na
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_getGeneHeatmapData
#'
.int_getGeneHeatmapData <- function(sfe, assay, subset_row, subset_col) {
  df <- SummarizedExperiment::assay(sfe, assay)

  if (!is.null(subset_row)) {
    df <- df[subset_row,]
  }
  if (!is.null(subset_col)) {
    df <- df[,subset_col]
  }

  return(df)
}


#' Internal: Get Annotations and Colours
#'
#' This function returns annotations label and colours for the heatmap
#'
#' @param fgwc gwc class object as generated by `fgwcSTE`
#' @param loc-annot which annotations to add
#' @param order_rows order rows based on annotation or cluster
#'
#' @returns a list
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_getAnnotsAndColours
#'
.int_getAnnotsAndColoursSFE <- function(sfe, loc_annot, order_col) {
  ## Check if loc_annot is empty
  if (isEmpty(loc_annot)) {
    if ("annotation" %in% colnames(colData(sfe))) {
      annotations <- colData(sfe)["annotation"]
    } else {
      annotations <- list()
    }
  } else {
    annotations <- as.data.frame(loc_annot)
    rownames(annotations) <- rownames(colData(sfe))
  }

  if (!is.null(order_col)) {
    annotations <- dplyr::arrange(annotations, .data[[order_col]])
  }

  annot_colours <- lapply(annotations, .int_fetchAnnotColours)

  list(annotations = annotations, colours = annot_colours)
}


#' Internal: Get Colours
#'
#' This function returns colours for the annotations. It is intended to work in
#' an `lapply` setup.
#'
#' @param annot annotations vector
#'
#' @returns a named vector of colours with annotation labels as names.
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_fetchAnnotColours
#'
.int_fetchAnnotColours <- function(annot) {
  ## A list of colours
  colour_list <- c("#1F78C8", "#FF0000", "#33A02C", "#6A33C2", "#FF7F00",
                   "#565656", "#FFD900", "#A6CEE3", "#FB6496", "#B2DF8A",
                   "#CAB2D6", "#FDBF6F", "#999999", "#EEE685", "#C8308C",
                   "#FFFFFF", "#720055", "#0000FF", "#36648B", "#00E2E5",
                   "#8B3B00", "#A52A3C", "#FF83FA", "#0000FF", "#00FF00",
                   "#000033", "#201A01", "#005300", "#FFC000", "#009FFF",
                   "#00FFBE", "#1F9698", "#B1CC71", "#F1085C", "#FE8F42",
                   "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F")
  ## Fetch unique annotation labels
  unique_annots <- unique(annot)
  ## Calculate the number of colours required
  colour_number <- length(unique_annots)
  ## Get colours
  annotation_colours <- colour_list[1:colour_number]
  ## Match colours to annotation labels
  names(annotation_colours) <- unique_annots
  ## Return named character string
  return(annotation_colours)
}






