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
#'                  analysis.
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
#' @importFrom purr reduce
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
                               type = c("spot", "hex"),
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
  }

  ## Fetch image if needed
  if (!is.null(res)) {
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
  ggplot(data) +
    tidyterra::geom_spatraster_rgb(data = image[[1]]) +
    ggplot2::geom_sf(data = data,
                     aes(geometry = geometry,
                         fill = expression),
                     alpha = alpha) +
    do.call(scale_fill_viridis_c, c(list(), fill_args)) +
    ggplot2::labs(fill = fill) +
    ggplot2::coord_sf() +
    ggplot2::facet_wrap(~gene, scales = "fixed", ncol = 1) +
    ggplot2::theme_void()

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
#'
.int_prepareDataMap <- function(sfe, assay, genes, minmax, type) {
  # genes_df <- assay(sfe, assay) %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   dplyr::select(tidyr::all_of(genes)) %>%
  #   tibble::rownames_to_column(var = "rownames")
  genes_df <- assay(sfe, assay) %>%
    t() %>%
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
