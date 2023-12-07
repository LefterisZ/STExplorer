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
  # stopifnot(is(sfe, "SpatialFeatureExperiment"))
  res <- match.arg(res)

  ## Select samples
  ids <- .int_getMSFEsmplID(list = msfe, sample_id = sample_id)

  ## Prepare neighbour graph data
  nb_data_list <-  lapply(msfe[ids], .int_getNBdata)
  ## Fetch image data and transform to raster
  if (plotImage) {
    image_list <- lapply(msfe[ids], .int_getImgDtMSFE, image_id = res)
  }
  ## Get capture area limits
  limits_list <- lapply(msfe[ids], .int_getImgLims)

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
#' @aliases .int_getImgDtMSFE
#' @importFrom SpatialFeatureExperiment imgRaster
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
    theme(plot.subtitle = element_text(hjust = 0.5))

  # return(p)
}
