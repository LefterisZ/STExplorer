#' Plot Quality Control of Spots in a SpatialFeatureExperiment
#'
#' This function plots the quality control of spots in a
#' SpatialFeatureExperiment, providing options for different spot types and
#' customization.
#'
#' @param m_sfe A SpatialFeatureExperiment or MetaSFE object.
#' @param type Character vector specifying the spot type, either "spot", "cntd"
#' or "hex".
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#' identifier(s). TRUE  equivalent to all samples/images, and NULL specifies
#' the first available entry.
#' @param in_tissue Logical, indicating whether to consider only spots in
#' tissue. Default is TRUE.
#' @param colours Vector of colours for spots. Default is NULL, which uses the
#' default colours.
#' @param y_reverse Logical, indicating whether to reverse the y-axis. Default
#' is FALSE.
#'
#' @return A ggplot object representing the quality control plot.
#'
#' @details
#' This function plots the quality control of spots in a
#' SpatialFeatureExperiment. It provides options for different spot types
#' (spot or hexagon or centroid) and customisation. The in_tissue argument
#' controls whether to consider only spots in tissue. The colours argument
#' allows users to customise spot colours, and y_reverse can be used to reverse
#' the y-axis.
#'
#' @seealso
#' Explore the related internal functions: \code{\link{.int_getSmplIDs}},
#' \code{\link{.int_plotSpot}}, and \code{\link{.int_plotHex}}.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial quality-control plotting ggplot2
#'
#' @rdname plotQC_spots
#' @aliases plotQC_spots
#'
#' @examples
#' # Example usage:
#' data(sfe)
#' plotQC_spots(sfe, type = "spot", sample_id = TRUE, in_tissue = TRUE,
#'              colours = NULL, y_reverse = FALSE)
#'
#' @export
plotQC_spots <- function(m_sfe,
                         type = c("spot", "hex", "cntd"),
                         sample_id = NULL,
                         in_tissue = TRUE,
                         colours = NULL,
                         y_reverse = FALSE){

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check arguments
  stopifnot(is.logical(in_tissue)) # check if logical is provided

  ## Check valid type argument
  type <- match.arg(type)
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
  }
  if (type == "spot") {
    stopifnot("spotPoly" %in% names(colGeometries(sfe)))
  }
  if (type == "cntd") {
    stopifnot("spotCntd" %in% names(colGeometries(sfe)))
  }

  ## Get the number of unique samples
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)
  n_samples <- length(unique(ids))

  ## Plot
  if (type == "spot") {
    p <- .int_plotSpot(sfe = sfe, in_tissue = in_tissue,
                       n_samples = n_samples, colours = colours, ids = ids)
  } else if (type == "hex") {
    p <- .int_plotHex(sfe = sfe, in_tissue = in_tissue,
                      n_samples = n_samples, colours = colours, ids = ids)
  } else if (type == "cntd") {
    p <- .int_plotCntd(sfe = sfe, in_tissue = in_tissue,
                       n_samples = n_samples, colours = colours, ids = ids)
  }

  ## Add facet_wrap if there are multiple samples
  if (n_samples > 1) {
    p <- p + ggplot2::facet_wrap(~ sample_id)
  }

  ## Remove strip from multifaceted plots
  if (in_tissue && n_samples > 1) {
    p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())
  }

  ## Reverse the y-axis if 'y_reverse' is TRUE
  if (y_reverse) {
    p <- p + ggplot2::scale_y_reverse()
  }

  ## Return the ggplot object
  return(p)

}


#' Plot Quality Control of Spots with Annotation in a SpatialFeatureExperiment
#'
#' This function plots the quality control of spots in a
#' SpatialFeatureExperiment, including annotation information.
#'
#' @param m_sfe A SpatialFeatureExperiment or MetaSFE object.
#' @param type Character vector specifying the spot type, either "spot", "cntd
#' or "hex".
#' @param fill_args List of arguments to customize the fill scale. Default is
#' an empty list.
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#' identifier(s). TRUE is equivalent to all samples/images, and NULL specifies
#' the first available entry.
#' @param colours Vector of colours for annotation. Default is NULL, which uses
#' default colours.
#' @param ... Additional arguments to be passed to the underlying plotting
#' functions.
#'
#' @return A ggplot object representing the quality control plot with
#' annotation.
#'
#' @details
#' This function plots the quality control of spots in a
#' SpatialFeatureExperiment, including annotation information.
#'
#' Users can specify the spot type (spots, centroids, or hexagons) and
#' customise the fill scale for annotation.
#'
#' The colours argument allows users to customize annotation colours, and
#' additional arguments (...) can be passed to the underlying plotting
#' functions.
#'
#' @seealso
#' Explore the related internal functions: \code{\link{.int_getSmplIDs}},
#' \code{\link{.int_dataToPlotAnnot}}, and \code{\link{getColours}}.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial quality-control annotation plotting ggplot2
#'
#' @rdname plotQC_spotsAnnotation
#' @aliases plotQC_spotsAnnotation
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(sfe)
#' plotQC_spotsAnnotation(sfe, type = "spot", fill_args = list(),
#'                        sample_id = TRUE, colours = NULL)
#' }
#'
#' @export
plotQC_spotsAnnotation <- function(m_sfe,
                                   sample_id = NULL,
                                   type = c("spot", "hex", "cntd"),
                                   fill_args = list(),
                                   colours = NULL,
                                   ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check valid type argument
  type <- match.arg(type)
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
  }
  if (type == "spot") {
    stopifnot("spotPoly" %in% names(colGeometries(sfe)))
  }
  if (type == "cntd") {
    stopifnot("spotCntd" %in% names(colGeometries(sfe)))
  }

  ## Get the number of unique samples
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)
  n_samples <- length(unique(ids))

  ## Fetch colours
  if (is.null(colours)) {
    annot_number <- length(unique(colData(sfe)$annotation))
    annot_cols <- getColours(annot_number)
  } else {
    annot_cols <- colours
  }

  ## Set some defaults if not provided
  if (isEmpty(fill_args)) {
    ## Make list
    fill_args <- list(values = annot_cols,
                      na.value = "grey95")
  }

  ## Put together the data to plot
  data <- .int_dataToPlotAnnot(sfe = sfe, type = type, ids = ids,
                               annotate = TRUE)

  ## Plot annotation map
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = data,
                     aes(geometry = geometry,
                         fill = annotation)) +
    do.call(scale_fill_manual, c(list(...), fill_args)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(fill = "Annotation")

  ## Add facet_wrap if there are multiple samples
  ## Or add plot title if there is only one plot
  if (n_samples > 1) {
    p <- p + ggplot2::facet_wrap(~ sample_id)
  } else {
    p <- p + ggplot2::labs(title = ids) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))
  }

  return(p)

}


#' Plot Tissue Image in a SpatialFeatureExperiment
#'
#' This function plots the tissue image in a SpatialFeatureExperiment,
#' providing options for different resolutions and spot types.
#'
#' @param m_sfe A SpatialFeatureExperiment or MetaSFE object.
#' @param res Character vector specifying the image resolution, either
#' "lowres", "hires", or "fullres".
#' @param type Character vector specifying the spot type, either "spot",
#' "hex", or "none".
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#' identifier(s). TRUE is equivalent to all samples/images, and NULL specifies
#' the first available entry.
#' @param fill_args List of arguments to customize the fill scale. Default is
#' an empty list.
#' @param annotate Logical, indicating whether to annotate the plot. Default
#' is FALSE.
#' @param colours Vector of colours for annotation. Default is NULL, which uses
#' default colours.
#' @param alpha Numeric, indicating the alpha (transparency) value for the
#' spots and the annotation (if selected). Default is 0.3.
#' @param ... Additional arguments to be passed to the underlying plotting functions.
#'
#' @details
#' This function plots the tissue image in a SpatialFeatureExperiment. Users
#' can specify the image resolution, spot type, and customize the fill scale
#' for annotation.
#' The annotate argument controls whether to annotate the plot, and colours
#' can be used to customise annotation colours. The alpha argument adjusts the
#' transparency of the spots overlaying the image.
#'
#' @seealso
#' Explore the related internal functions: \code{\link{.int_getImgDt}},
#' \code{\link{.int_getImgLims}}, and \code{\link{.int_tissueImgPlot}}.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial quality-control annotation plotting ggplot2
#'
#' @importFrom patchwork wrap_plots
#'
#' @rdname plotQC_tissueImg
#' @aliases plotQC_tissueImg
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(sfe)
#' # With annotations
#' plotQC_tissueImg(sfe, res = "lowres", type = "spot", annotate = TRUE)
#'
#' # Without annotations
#' plotQC_tissueImg(sfe, res = "lowres", type = "spot", annotate = FALSE)
#' }
#'
#' @export
plotQC_tissueImg <- function(m_sfe,
                             res = c("lowres", "hires", "fullres"),
                             type = c("spot", "hex", "cntd", "none"),
                             sample_id = NULL,
                             fill_args = list(),
                             annotate = FALSE,
                             colours = NULL,
                             alpha = 0.3,
                             ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check valid type argument
  type <- match.arg(type)
  res <- match.arg(res)
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
  }
  if (type == "spot") {
    stopifnot("spotPoly" %in% names(colGeometries(sfe)))
  }
  if (type == "cntd") {
    stopifnot("spotCntd" %in% names(colGeometries(sfe)))
  }

  ## Fetch image data and transform to raster
  image_list <- .int_getImgDt(sfe = sfe, sample_id = sample_id, image_id = res)
  ## Get capture area limits
  limits_list <- .int_getImgLims(sfe = sfe)

  ## Plot using ggplot2
  plots <- lapply(seq_along(image_list), function(i) {
    .int_tissueImgPlot(image = image_list[[i]],
                       image_name = names(image_list)[i],
                       limits_list = limits_list,
                       sfe = sfe,
                       type = type,
                       annotate = annotate,
                       colours = colours,
                       fill_args = fill_args,
                       alpha = alpha,
                       ...)
  })

  ## Combine the individual plots into a patchwork
  patchwork::wrap_plots(plots)

  # return(p)
}


#' Plot Quality Control Histogram in a SpatialFeatureExperiment
#'
#' This function plots the quality control histogram in a
#' SpatialFeatureExperiment, providing options for different metrics,
#' sample IDs, and customisation.
#'
#' @param m_sfe A SpatialFeatureExperiment or MetaSFE object.
#' @param limits Numeric vector specifying the limits of the x-axis.
#' @param metric Character vector specifying the metric, either "libsize",
#' "mito", "detected", "sizeFact", or "cellCount".
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#' identifier(s). TRUE is equivalent to all samples/images, and NULL specifies
#' the first available entry.
#' @param hist_args List of arguments to customise the histogram. Default is an
#' empty list.
#' @param dens_args List of arguments to customise the density plot. Default is
#' an empty list.
#' @param vline_args List of arguments to customise the vertical line. Default
#' is an empty list.
#' @param theme_args List of arguments to customise theme elements.
#' @param ... Additional arguments to be passed to the underlying plotting
#' functions.
#'
#' @return A ggplot object representing the quality control histogram.
#'
#' @details
#' This function plots the quality control histogram in a
#' SpatialFeatureExperiment. Users can specify the metric (libsize, mito,
#' detected, sizeFact, or cellCount), sample ID, and customize the appearance
#' of the histogram, density plot, and vertical line.
#'
#' @seealso
#' Explore the related internal functions: \code{\link{.int_getSmplIDs}}.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial quality-control histogram plotting ggplot2
#'
#' @importFrom ggplot2 geom_histogram geom_density geom_vline after_stat
#'
#' @rdname plotQC_hist
#' @aliases plotQC_hist
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(sfe)
#' plotQC_hist(sfe, limits = c(0, 5000), metric = "libsize")
#' }
#'
#' @export
plotQC_hist <- function(m_sfe,
                        limits,
                        metric = c("libsize", "mito",
                                   "detected", "cellCount"),
                        sample_id = TRUE,
                        hist_args = list(),
                        dens_args = list(),
                        vline_args = list(),
                        theme_args = list(),
                        ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get the number of unique samples
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)
  n_samples <- length(unique(ids))

  ## Check valid type argument
  metric <- match.arg(metric)

  ## Fetch the data
  if (metric == "libsize") {
    aes_x <- sfe$sum
    lab_x <- "Library Size"
  } else if (metric == "mito") {
    aes_x <- sfe$subsets_mito_percent
    lab_x <- "Percentage of mitochondrial expression"
  } else if (metric == "detected") {
    aes_x <- sfe$detected
    lab_x <- "Genes expressed in each spot"
  } else if (metric == "sizeFact") {
    aes_x <- sfe$sizeFactor
    lab_x <- "Size Factor"
  } else if (metric == "cellCount") {
    aes_x <- sfe$cellCount
    lab_x <- "Number of Cells"
  }
  data <- data.frame(aes_x = aes_x,
                     sample_id = sfe$sample_id) %>%
    dplyr::filter(.data$sample_id %in% ids)

  ## Set some defaults if not provided
  if (DelayedArray::isEmpty(hist_args)) {
    hist_args <- list(colour = "black",
                      fill = "grey",
                      bins = 50)
  }
  if (DelayedArray::isEmpty(dens_args)) {
    dens_args <- list(alpha = 0.5,
                      adjust = 0.5,
                      fill = "#A0CBE8",
                      colour = "#4E79A7")
  }
  if (DelayedArray::isEmpty(vline_args)) {
    vline_args <- list(colour = "red",
                       linetype = "dashed")
  }
  if (missing(limits)) {
    limits = c(0, max(aes_x))
  }

  ## Generate plot
  p <- ggplot(data = data,
         aes(x = aes_x)) +
    do.call(geom_histogram, c(list(...), hist_args,
                              list(aes(y = after_stat(density))))) +
    do.call(geom_density, c(list(...), dens_args)) +
    do.call(geom_vline, c(list(...), vline_args, list(xintercept = limits))) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::xlab(lab_x) +
    ggplot2::ylab("Density") +
    ggplot2::theme_classic()

  ## Add facet_wrap if there are multiple samples
  if (n_samples > 1) {
    p <- p + ggplot2::facet_wrap(~ sample_id)
  }

  ## Modify theme according to user input (if present)
  if (DelayedArray::isEmpty(theme_args)) {
    p <- p + do.call(theme, c(hist_args))
  }

  return(p)
}


#' Plot Quality Control Scatter Plot in a SpatialFeatureExperiment
#'
#' This function plots the quality control scatter plot in a
#' SpatialFeatureExperiment, providing options for different metrics,
#' sample IDs, and customisation.
#'
#' @param m_sfe A SpatialFeatureExperiment or MetaSFE object.
#' @param limits Numeric vector specifying the limits of the x and y axes.
#' @param metric Character vector specifying the metric, either "libsize",
#' "mito", or "detected".
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#' identifier(s). TRUE is equivalent to all samples/images, and NULL specifies
#' the first available entry.
#' @param point_args List of arguments to customize the scatter points.
#' Default is an empty list.
#' @param hline_args List of arguments to customize the horizontal line.
#' Default is an empty list.
#' @param ... Additional arguments to be passed to the underlying plotting
#' functions.
#'
#' @return A ggplot object representing the quality control scatter plot.
#'
#' @details
#' This function plots the quality control scatter plot in a
#' SpatialFeatureExperiment. Users can specify the metric (libsize, mito, or
#' detected), sample ID, and customize the appearance of the scatter points
#' and horizontal line.
#'
#' @seealso
#' Explore the related internal functions: \code{\link{.int_getSmplIDs}}.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial quality-control scatter-plot plotting ggplot2
#'
#' @rdname plotQC_scat
#' @aliases plotQC_scat
#'
#' @importFrom ggExtra ggMarginal
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(sfe)
#' plotQC_scat(sfe, limits = c(0, 5000), metric = "libsize")
#' }
#'
#' @export
plotQC_scat <- function(m_sfe,
                        limits,
                        metric = c("libsize", "mito",
                                   "detected", "cellCount"),
                        sample_id = TRUE,
                        point_args = list(),
                        hline_args = list(),
                        ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check valid type argument
  metric <- match.arg(metric)

  ## Get the number of unique samples
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)
  n_samples <- length(unique(ids))

  ## Fetch the data
  if (metric == "libsize") {
    aes_y <- sfe$sum
    lab_y <- "Library Size"
  } else if (metric == "mito") {
    aes_y <- sfe$subsets_mito_percent
    lab_y <- "Percentage of mitochondrial expression"
  } else if (metric == "detected") {
    aes_y <- sfe$detected
    lab_y <- "Genes expressed in each spot"
  }
  data <- data.frame(x = sfe$cellCount,
                     y = aes_y,
                     sample_id = sfe$sample_id) %>%
    dplyr::filter(sample_id %in% ids)

  ## Set some defaults if not provided
  if (DelayedArray::isEmpty(point_args)) {
    point_args <- list(colour = "black",
                      size = 2)
  }
  if (DelayedArray::isEmpty(hline_args)) {
    hline_args <- list(colour = "red",
                       linetype = "dashed")
  }
  if (missing(limits)) {
    limits = c(0, max(aes_y))
  }

  ## Generate plot
  p <- ggplot(data = data,
         aes(x = data$x, y = data$y)) +
    do.call(geom_point, c(list(...), point_args)) +
    do.call(geom_hline, c(list(...), hline_args, list(yintercept = limits))) +
    ggplot2::geom_smooth(se = FALSE) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::xlab("Number of Cells") +
    ggplot2::ylab(lab_y) +
    ggplot2::theme_classic()

  p <- ggExtra::ggMarginal(p, type = 'histogram', margins = 'both')

  return(p)
}


#' Plot Quality Control Map in a SpatialFeatureExperiment
#'
#' This function plots the quality control map plot in a
#' SpatialFeatureExperiment, providing options for different metrics,
#' sample IDs, and customisation.
#'
#' @param m_sfe An object of class 'SpatialFeatureExperiment' or
#'              'MetaSpatialFeatureExperiment' containing spatial
#'              data.
#' @param sample_id A character string specifying the sample ID for
#'                  analysis. If NULL the first sample is selected.
#' @param type Character vector specifying the spot type, either "spot", "cntd"
#'             or "hex".
#' @param metric Character vector specifying the metric, either "libsize",
#'               "mito", "detected", or custom. If custom is selected then you
#'               have to provide the `metric_name` and `metric_lab` arguments
#'               too.
#' @param colours One of "viridis" or "custom". Select between Viridis' package
#'                colour schemes and custom colours provided by you. IMPORTANT:
#'                Look at Details for more information.
#' @param col_args List of arguments to customise the map fill colour.
#'                 Default is an empty list. If you select "viridis" at the
#'                 `colours` argument, here you can select the colour scheme
#'                 you prefer. I.e., `list(option = "magma")`.
#' @param metric_name Character string. The custom metric's column name as
#'                    found in the `colData`.
#' @param metric_lab Character string. The title you want for the map's legend.
#' @param ... Additional arguments to be passed to the underlying `geom_sf`
#'            plotting function.
#'
#' @return A ggplot object representing the quality control map plot.
#'
#' @details
#' This function plots the quality control map plot from a
#' SpatialFeatureExperiment. Users can specify the metric (libsize, mito,
#' cellCount, detected, or custom), sample ID, and customise the appearance of
#' the map's points.
#'
#' About using custom colours:
#' At the moment, the function is implementing `scale_fill_gradientn` which
#' expects a gradient like the ones generated by `colorRampPalette`. For
#' example: `colors <- colorRampPalette(c("blue", "red"))(100)`
#'
#' @seealso
#' Explore the related internal functions: \code{\link{.int_getSmplIDs}}.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial quality-control scatter-plot plotting ggplot2
#'
#' @rdname plotQC_map
#' @aliases plotQC_map
#'
#' @importFrom DelayedArray isEmpty
#' @importFrom ggplot2 scale_colour_viridis_c scale_fill_viridis_c
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(sfe)
#' plotQC_map(sfe, metric = "libsize")
#' }
#'
#' @export
plotQC_map <- function(m_sfe,
                       metric = c("libsize", "mito",
                                  "detected", "cellCount", "custom"),
                       sample_id = NULL,
                       type = c("spot", "hex", "cntd"),
                       colours = c("viridis", "custom"),
                       col_args = list(),
                       metric_name = NULL,
                       metric_lab = NULL,
                       ...) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check type argument
  if (missing(type)) {
    type <- "hex"
  }
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
    type <- "spotHex"
  }
  if (type == "spot") {
    stopifnot("spotPoly" %in% names(colGeometries(sfe)))
    type <- "spotPoly"
  }
  if (type == "cntd") {
    stopifnot("spotCntd" %in% names(colGeometries(sfe)))
    type <- "spotCntd"
  }
  ## Check valid metric argument
  metric <- match.arg(metric)

  ## Fetch the data
  if (metric == "libsize") {
    fill <- sfe$sum
    lab_fill <- "Library Size"
  } else if (metric == "mito") {
    fill <- sfe$subsets_mito_percent
    lab_fill <- "% mitochondrial\nexpression"
  } else if (metric == "detected") {
    fill <- sfe$detected
    lab_fill <- "No. of Genes"
  } else if (metric == "cellCount") {
    fill <- sfe$cell_count
    lab_fill <- "No. of Cells"
  } else if (metric == "custom") {
    fill <- sfe[[metric_name]]
    lab_fill <- metric_lab
  }
  data <- data.frame(fill = fill,
                     geometry = colGeometry(sfe, type = type)[["geometry"]])

  ## Set some defaults if not provided
  if (missing(colours)) {
    colours <- "viridis"
  }
  if (DelayedArray::isEmpty(col_args) & colours == "viridis") {
    col_args <- list(option = "viridis")
  } else if (DelayedArray::isEmpty(col_args) & colours == "custom") {
    stop("When selecting custom colours, `col_args` must not be empty.")
  }

  ## Generate plot
  if (type %in% c("spotPoly", "spotHex")) {
    p <- ggplot(data) + ggplot2::geom_sf(aes(geometry = geometry,
                                             fill = fill),
                                         ...) +
      ggplot2::labs(fill = lab_fill)
  } else {
    p <- ggplot(data) + ggplot2::geom_sf(aes(geometry = geometry,
                                             colour = fill),
                                         ...) +
      ggplot2::labs(colour = lab_fill)
  }

  if (colours == "viridis") {
    if (type == "spotCntd") {
      p <- p + do.call(scale_colour_viridis_c, c(list(), col_args))
    } else {
      p <- p + do.call(scale_fill_viridis_c, c(list(), col_args))
    }

  } else if (colours == "custom") {
    if (type == "spotCntd") {
      p <- p + do.call(scale_colour_gradientn, c(list(), col_args))
    } else {
      p <- p + do.call(scale_fill_gradientn, c(list(), col_args))
    }
  }

   p +
     ggplot2::coord_sf() +
     ggplot2::theme_void()
}


#' Plot Quality Control Filtered Locations Map in a SpatialFeatureExperiment
#'
#' This function plots the quality control filtered locations map in a
#' SpatialFeatureExperiment, providing options for different metrics and
#' sample IDs.
#'
#' @param m_sfe A SpatialFeatureExperiment or MetaSFE object.
#' @param metric Character vector specifying the metric, either "libsize",
#' "mito", "NAs", "detected", "cellCount", "discard", or custom. If custom is
#' selected then you have to provide the `metric_name` and `metric_lab`
#' arguments too.
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#' identifier(s). TRUE is equivalent to all samples/images, and NULL specifies
#' the first available entry.
#' @param type Character vector specifying the spot type, either "spot", "cntd"
#'             or "hex".
#' @param metric_name Character string. The custom metric's column name as
#'                    found in the `colData`.
#' @param metric_lab Character string. The title you want for the map's legend.
#' @param ... Additional arguments to be passed to the underlying `geom_sf`
#'            plotting function.
#'
#' @details
#' This function plots the quality control filtered locations map in a
#' SpatialFeatureExperiment. Users can specify the metric (libsize, mito, NAs,
#' detected, cellCount, or discard) and sample ID.
#'
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial quality-control filtering mapping ggplot2
#'
#' @rdname plotQC_filtered
#' @aliases plotQC_filtered
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(sfe)
#' plotQC_filtered(sfe, metric = "libsize")
#' }
#'
#' @export
plotQC_filtered <- function(m_sfe,
                            metric = c("libsize", "mito", "NAs", "detected",
                                       "cellCount", "discard", "custom"),
                            sample_id = TRUE,
                            type = c("spot", "hex", "cntd"),
                            metric_name = NULL,
                            metric_lab = NULL,
                            ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check type argument
  if (missing(type)) {
    type <- "hex"
  }
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
    type <- "spotHex"
  }
  if (type == "spot") {
    stopifnot("spotPoly" %in% names(colGeometries(sfe)))
    type <- "spotPoly"
  }
  if (type == "cntd") {
    stopifnot("spotCntd" %in% names(colGeometries(sfe)))
    type <- "spotCntd"
  }

  ## Check valid metric argument
  metric <- match.arg(metric)

  ## Get the number of unique samples
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)
  n_samples <- length(unique(ids))

  ## Prepare data for plotting
  if (metric == "libsize") {
    fill <- "qc_lib_size"
    title <- "Filtered for Library Size"
  } else if (metric == "mito") {
    fill <- "qc_mito"
    title <- "Filtered for % of Mitochondrial Expression"
  } else if (metric == "detected") {
    fill <- "qc_detected"
    title <- "Filter for Number of Genes"
  } else if (metric == "cellCount") {
    fill <- "qc_cellCount"
    title <- "Filter for number of cells"
  } else if (metric == "NAs") {
    fill <- "qc_NA_spots"
    title <- "Filter for locations without annotation"
  } else if (metric == "discard") {
    fill <- "qc_discard"
    title <- "Overall Low Quality"
  } else if (metric == "custom") {
    fill <- metric_name
    title <- metric_lab
  }

  if (metric == "discard") {
    fill_label <- "Discarded"
  } else {
    fill_label <- "To be \nDiscarded"
  }

  ## Put together the data to plot
  data <- .int_dataToPlot(sfe = sfe, ids = ids, fill = fill, type = type)

  ## Plot filtered locations map
  p <- ggplot2::ggplot()

  if (type == "spotCntd") {
    p <- p +
      ggplot2::geom_sf(data = data,
                       ggplot2::aes(geometry = geometry,
                                    colour = fill)) +
      ggplot2::scale_colour_manual(values = c("grey95", "red")) +
      ggplot2::labs(colour = fill_label,
                    title = title)
  } else {
    p <- p +
      ggplot2::geom_sf(data = data,
                       ggplot2::aes(geometry = geometry,
                                    fill = fill)) +
      ggplot2::scale_fill_manual(values = c("grey95", "red")) +
      ggplot2::labs(fill = fill_label,
                    title = title)
  }

  p <- p +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")

  ## Add facet_wrap if there are multiple samples
  ## Or add plot title if there is only one plot
  if (n_samples > 1) {
    p <- p + ggplot2::facet_wrap(~ sample_id)
  } else {
    p <- p + ggplot2::labs(subtitle = ids) +
      ggplot2::theme(plot.subtitle = element_text(hjust = 0.5))
  }

  return(p)
}


#' Plot Size Factors Histogram and Density in a SpatialFeatureExperiment
#'
#' This function plots the size factors histogram and density in a
#' SpatialFeatureExperiment, providing options for sample IDs and customisation.
#'
#' @param m_sfe A SpatialFeatureExperiment or MetaSFE object.
#' @param sample_id Character string, TRUE, or NULL specifying sample/image
#' identifier(s). TRUE is equivalent to all samples/images, and NULL specifies
#' the first available entry.
#' @param hist_args List of arguments to customize the histogram. Default is
#' an empty list.
#' @param dens_args List of arguments to customize the density plot. Default is
#' an empty list.
#' @param ... Additional arguments to be passed to the underlying plotting
#' functions.
#'
#' @return A ggplot object representing the size factors histogram and density.
#'
#' @details
#' This function plots the size factors histogram and density in a
#' SpatialFeatureExperiment. Users can specify the sample ID and customise the
#' appearance of the histogram and density plot.
#'
#' @seealso
#' Explore the related internal functions: \code{\link{.int_getSmplIDs}}.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial size-factors histogram density plotting ggplot2
#'
#' @rdname plotQC_sizeFactors
#' @aliases plotQC_sizeFactors
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(sfe)
#' plotQC_sizeFactors(sfe)
#' }
#'
#' @export
plotQC_sizeFactors <- function(m_sfe,
                               sample_id = TRUE,
                               hist_args = list(),
                               dens_args = list(),
                               ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Get the number of unique samples (in case it is an SFE with multiple samples)
  ids <- .int_getSmplIDs(sfe = sfe, sample_id = sample_id)
  n_samples <- length(unique(ids))

  ## Fetch the data
  data <- data.frame(aes_x = sfe$sizeFactor,
                     sample_id = sfe$sample_id) %>%
    dplyr::filter(.data$sample_id %in% ids)
  lab_x <- "Size Factor"

  ## Set some defaults if not provided
  if (DelayedArray::isEmpty(hist_args)) {
    hist_args <- list(colour = "black",
                      fill = "grey",
                      bins = 50)
  }
  if (DelayedArray::isEmpty(dens_args)) {
    dens_args <- list(alpha = 0.5,
                      adjust = 0.5,
                      fill = "#A0CBE8",
                      colour = "#4E79A7")
  }

  ## Generate plot
  p <- ggplot(data = data,
         aes(x = aes_x)) +
    do.call(geom_histogram, c(list(...), hist_args,
                              list(aes(y = after_stat(density))))) +
    do.call(geom_density, c(list(...), dens_args)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::xlab(lab_x) +
    ggplot2::ylab("Density") +
    ggplot2::theme_classic()

  ## Add facet_wrap if there are multiple samples
  if (n_samples > 1) {
    p <- p + ggplot2::facet_wrap(~ sample_id)
  }

  return(p)
}


# ---------------------------------------------------------------------------- #
#  ############### INTERNAL FUNCTIONS ASSOCIATED WITH QC PLOTS ###############
# ---------------------------------------------------------------------------- #
#' Internal Function: .int_plotSpot
#'
#' Generates a ggplot2 plot of spatial coordinates for spots with
#' optional colour customisation.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param in_tissue Logical. If TRUE, plot only in-tissue spots; if FALSE,
#' plot all spots.
#' @param n_samples Number of samples in the dataset.
#' @param colours Vector of colours for filling the spots. If NULL, default
#' colours are used.
#' @param ids Vector of sample IDs to include in the plot.
#'
#' @return Returns a ggplot2 plot object.
#'
#' @details This function combines sample metadata with spatial coordinates and
#' creates a plot using ggplot2.
#' The plot can be customised based on the input parameters.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_plotSpot
#'
.int_plotSpot <- function(sfe, in_tissue, n_samples, colours, ids) {
  ## Put together the data
  data <- data.frame(sample_id = colData(sfe)$sample_id,
                     in_tissue = colData(sfe)$in_tissue,
                     geometry = colGeometries(sfe)$spotPoly$geometry) %>%
    dplyr::filter(.data$sample_id %in% ids)

  ## Select arguments and data for in-tissue spots plot only or all spots
  if (in_tissue) {
    data <- data[data$in_tissue, ]
    fill <- "Sample ID"
    geomSF_args <- list(aes(geometry = data$geometry, fill = data$sample_id))
    if (is.null(colours)) {
      colours <- getColours(n_samples)
    } else {
      colours <- colours
    }

  } else {
    fill <- "In tissue"
    geomSF_args <- list(aes(geometry = data$geometry, fill = data$in_tissue))
    if (is.null(colours)) {
      colours <- c("#1F78C8", "#FF0000")
    } else {
      colours <- colours
    }
  }

  ## Create the ggplot object

  ## To silence the 'coordinate system already present' warning as per:
  ## https://github.com/tidyverse/ggplot2/issues/2799
  cf <- coord_fixed()
  cf$default <- TRUE

  p <- ggplot2::ggplot() +
    ggplot2::scale_fill_manual(values = colours) +
    cf +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = "Spatial coordinates",
                  fill = fill)

  ## Add geom_sf
  ## - for an unknown reason, when the geom_sf() is added inside the above it
  ##   throws an error
  p <- p + do.call(geom_sf, c(list(data = data, lwd = 0), geomSF_args))

  return(p)
}


#' Internal Function: .int_plotCntd
#'
#' Generates a ggplot2 plot of spatial coordinates for spots with
#' optional colour customisation.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param in_tissue Logical. If TRUE, plot only in-tissue spots; if FALSE,
#' plot all spots.
#' @param n_samples Number of samples in the dataset.
#' @param colours Vector of colours for filling the spots. If NULL, default
#' colours are used.
#' @param ids Vector of sample IDs to include in the plot.
#'
#' @return Returns a ggplot2 plot object.
#'
#' @details This function combines sample metadata with spatial coordinates and
#' creates a plot using ggplot2.
#' The plot can be customised based on the input parameters.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_plotCntd
#'
.int_plotCntd <- function(sfe, in_tissue, n_samples, colours, ids) {
  ## Put together the data
  data <- data.frame(sample_id = colData(sfe)$sample_id,
                     in_tissue = colData(sfe)$in_tissue,
                     geometry = colGeometries(sfe)$spotCntd$geometry) %>%
    dplyr::filter(.data$sample_id %in% ids)

  ## Select arguments and data for in-tissue spots plot only or all spots
  if (in_tissue) {
    data <- data[data$in_tissue, ]
    colour <- "Sample ID"
    geomSF_args <- list(aes(geometry = data$geometry, colour = data$sample_id))
    if (is.null(colours)) {
      colours <- getColours(n_samples)
    } else {
      colours <- colours
    }

  } else {
    colour <- "In tissue"
    geomSF_args <- list(aes(geometry = data$geometry, colour = data$in_tissue))
    if (is.null(colours)) {
      colours <- c("#1F78C8", "#FF0000")
    } else {
      colours <- colours
    }
  }

  ## Create the ggplot object

  ## To silence the 'coordinate system already present' warning as per:
  ## https://github.com/tidyverse/ggplot2/issues/2799
  cf <- coord_fixed()
  cf$default <- TRUE

  p <- ggplot2::ggplot() +
    ggplot2::scale_colour_manual(values = colours) +
    cf +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = "Spatial coordinates",
                  colour = colour)

  ## Add geom_sf
  ## - for an unknown reason, when the geom_sf() is added inside the above it
  ##   throws an error
  p <- p + do.call(geom_sf, c(list(data = data, lwd = 0), geomSF_args))

  return(p)
}


#' Internal Function: .int_plotHex
#'
#' Generates a ggplot2 plot of spatial coordinates for hexagonal spots with
#' optional colour customisation.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param in_tissue Logical. If TRUE, plot only in-tissue spots; if FALSE,
#' ]plot all spots.
#' @param n_samples Number of samples in the dataset.
#' @param colours Vector of colours for filling the hexagonal spots. If NULL,
#' default colours are used.
#' @param ids Vector of sample IDs to include in the plot.
#'
#' @return Returns a ggplot2 plot object.
#'
#' @details This function combines sample metadata with spatial coordinates
#' and creates a plot using ggplot2.
#' The plot can be customized based on the input parameters.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_plotHex
#'
.int_plotHex <- function(sfe, in_tissue, n_samples, colours, ids) {
  ## Combine sample metadata with spatial coordinates
  data <- data.frame(sample_id = colData(sfe)$sample_id,
                     in_tissue = colData(sfe)$in_tissue,
                     geometry = colGeometries(sfe)$spotHex$geometry) %>%
    dplyr::filter(.data$sample_id %in% ids)

  ## Select arguments and data for in-tissue spots plot only or all spots
  if (in_tissue) {
    data <- data[data$in_tissue, ]
    fill <- "Sample ID"
    geomSF_args <- list(aes(geometry = data$geometry, fill = data$sample_id))
    if (is.null(colours)) {
      colours <- getColours(n_samples)
    } else {
      colours <- colours
    }

  } else {
    fill <- "In tissue"
    geomSF_args <- list(aes(geometry = data$geometry, fill = data$in_tissue))
    if (is.null(colours)) {
      colours <- c("#1F78C8", "#FF0000")
    } else {
      colours <- colours
    }
  }

  ## Create the ggplot object

  ## To silence the 'coordinate system already present' warning as per:
  ## https://github.com/tidyverse/ggplot2/issues/2799
  cf <- coord_fixed()
  cf$default <- TRUE

  p <- ggplot2::ggplot() +
    ggplot2::scale_fill_manual(values = colours) +
    cf +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = "Spatial coordinates",
                  fill = fill)

  ## Add geom_sf
  ## - for an unknown reason, when the geom_sf() is added inside the above it
  ##   throws an error
  p <- p + do.call(geom_sf, c(list(data = data), geomSF_args))

  return(p)
}


#' Internal Function: .int_dataToPlot
#'
#' Creates a data frame for generating spatial plots based on specified fill
#' variable.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param ids Vector of sample IDs to include in the plot.
#' @param fill Variable used for filling the spots in the plot.
#' @param type Acharacter string indicating the geometries to use.
#'
#' @return Returns a data frame suitable for generating spatial plots.
#'
#' @details This function creates a data frame with sample metadata and spatial
#' coordinates based on the specified fill variable.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_dataToPlot
#'
#' @importFrom dplyr filter
#'
.int_dataToPlot <- function(sfe, ids, fill, type) {
  ## Create a data frame
  data <- data.frame(sample_id = colData(sfe)$sample_id,
                     fill = colData(sfe)[[fill]])

  ## Add geometries according to requested type
  data$geometry <- colGeometry(sfe, type = type)[["geometry"]]

  data <- data %>%
    dplyr::filter(.data$sample_id %in% ids)

  return(data)
}


#' Internal Function: .int_dataToPlotAnnot
#'
#' Creates a data frame for generating spatial plots with optional annotations.
#'
#' @param sfe An object containing spatial feature data.
#' @param type Character vector specifying the type of plot ("spot" or "hex").
#' @param ids Vector of sample IDs to include in the plot.
#' @param annotate Logical. If TRUE, includes annotations in the data frame.
#'
#' @return Returns a data frame suitable for generating spatial plots with
#' annotations.
#'
#' @details This function creates a data frame with sample metadata, spatial
#' coordinates, and optional annotations.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_dataToPlotAnnot
#'
#' @importFrom dplyr filter
#'
.int_dataToPlotAnnot <- function(sfe, type = c("spot", "hex"), ids, annotate) {
  ## Create a data frame
  if (annotate) {
    data <- data.frame(sample_id = colData(sfe)$sample_id,
                       in_tissue = colData(sfe)$in_tissue,
                       annotation = colData(sfe)$annotation)
  } else {
    data <- data.frame(sample_id = colData(sfe)$sample_id,
                       in_tissue = colData(sfe)$in_tissue)
  }


  ## Add geometries according to requested type
  data$geometry <- .int_selectGeom(sfe, type = type)
  data <- data[data$in_tissue, ] %>% # in case raw feature matrix was imputed
    dplyr::filter(.data$sample_id %in% ids)

  return(data)
}


#' Internal: Select geometry
#'
#' An internal function to select geometries from a SpatialFeatureExperiment
#' object based on the specified type.
#'
#' @author Eleftherios Zormpas
#' @param sfe An object of class SpatialFeatureExperiment.
#' @param type A character string specifying the type of geometry to select.
#' Options are "spot" (spot geometry), "hex" (hexagon geometry), or "centroid"
#' (centroid geometry). Default is "spot".
#'
#' @return A geometry object representing the selected geometry type.
#'
#' @details This function retrieves the specified type of geometry from a
#' SpatialFeatureExperiment object. The available options for the type
#' parameter are "spot" (spot geometry), "hex" (hexagon geometry), and
#' "centroid" (centroid geometry).
#'
.int_selectGeom <- function(sfe, type = c("spot", "hex", "centroid")) {
  ## Add geometries according to requested type
  if (type == "spot") {
    return(colGeometries(sfe)$spotPoly$geometry)
  } else if (type == "hex") {
    return(colGeometries(sfe)$spotHex$geometry)
  } else if (type == "centroid") {
    return(colGeometries(sfe)$spotCntd$geometry)
  }
}


#' Internal Function: .int_annotColours
#'
#' Fetches annotation colours for a SpatialFeatureExperiment object.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param colours Vector of colours for annotations. If NULL, default colours
#' are used.
#'
#' @return Returns a vector of annotation colours.
#'
#' @details This function retrieves annotation colours based on the provided
#' colours vector or generates default colours if none are provided.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_annotColours
#'
.int_annotColours <- function(sfe, colours) {
  ## Fetch colours
  if (is.null(colours)) {
    annot_number <- length(unique(colData(sfe)$annotation))
    annot_cols <- getColours(annot_number)
  } else {
    annot_cols <- colours
  }

  return(annot_cols)
}


#' Internal Function: .int_fillArgs
#'
#' Sets default fill arguments if not provided.
#'
#' @param fill_args List of fill arguments.
#' @param annot_cols Vector of annotation colours.
#'
#' @return Returns a list of fill arguments with defaults if not provided.
#'
#' @details This function sets default fill arguments, such as values and
#' na.value, if the input list is empty.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_fillArgs
#'
.int_fillArgs <- function(fill_args, annot_cols) {
  ## Set some defaults if not provided
  if (isEmpty(fill_args)) {
    ## Make list
    fill_args <- list(values = annot_cols,
                      na.value = "grey95")
  }

  return(fill_args)
}


#' Internal Function: .int_getImgLims
#'
#' Extracts image limits for plotting based on spot locations.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param res The resolution as provided in the main function.
#'
#' @return Returns a list of x and y limits for image plotting.
#'
#' @details This function extracts spot locations and sets x and y limits for image plotting.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_getImgLims
#'
.int_getImgLims <- function(sfe, res) {
  ## Extract the spot locations
  spot_coords <- spatialCoords(sfe) %>% as.data.frame()

  ## Set limits
  xlim <- c(min(spot_coords$pxl_col_in_fullres) - 100,
            max(spot_coords$pxl_col_in_fullres) + 100)
  ylim <- c(min(spot_coords$pxl_row_in_fullres) - 100,
            max(spot_coords$pxl_row_in_fullres) + 100)

  limits_list <- list()
  limits_list[[1]] <- xlim
  limits_list[[2]] <- ylim

  return(limits_list)
}


#' Internal Function: .int_getImgSmplIDs
#'
#' Fetches sample IDs for image plotting.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#'
#' @return Returns a character vector of sample IDs for image plotting.
#'
#' @details This function fetches sample IDs from image data based on the input
#' parameters, providing flexibility for customising the sample IDs for image
#' plotting.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_getImgSmplIDs
#'
#' @importFrom SpatialFeatureExperiment imgData
#'
.int_getImgSmplIDs <- function(sfe, sample_id = NULL) {
  ## Fetch the required sample IDs
  if (is.null(sample_id)) {
    sample_id <- imgData(sfe)$sample_id[1]
  } else if (is.character(sample_id)) {
    sample_id <- sample_id
  } else if (sample_id) {
    sample_id <- imgData(sfe)$sample_id
  }

  return(sample_id)
}


#' Internal Function: .int_getImgDt
#'
#' Fetches image data and transforms it to raster for plotting.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#' @param image_id Character vector specifying the image type ("lowres" by
#' default).
#'
#' @return Returns a list of raster objects based on image data.
#'
#' @details This function fetches image data, transforms it to raster, and
#' organises it into a list with sample IDs as names.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @importFrom SpatialFeatureExperiment getImg
#'
#' @rdname dot-int_getImgDt
#'
.int_getImgDt <- function(sfe, sample_id = NULL, image_id = "lowres") {
  ## Fetch image data and transform to raster
  image_list <- getImg(sfe, sample_id = sample_id, image_id = image_id)
  if (is.list(image_list)) {
    raster_list <- lapply(image_list, imgRaster)
  } else {
    raster_list <- list()
    raster_list[[1]] <- imgRaster(image_list)
  }

  ## Fetch the required sample IDs
  sample_id <- .int_getImgSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Add names to list
  names(raster_list) <- sample_id

  return(raster_list)
}


#' Internal Function: .int_tissueImgPlot
#'
#' Generates a ggplot2 plot for tissue images with optional annotations.
#'
#' @param image List of raster objects.
#' @param image_name Character vector specifying the image name.
#' @param limits_list List of x and y limits for image plotting.
#' @param sfe A SpatialFeatureExperiment object.
#' @param type Character vector specifying the plot type ("none", "spot", or
#' "hex").
#' @param annotate Logical. If TRUE, includes annotations in the plot.
#' @param colours Vector of colours for annotations.
#' @param fill_args List of fill arguments.
#' @param alpha Numeric. Transparency level for image plotting.
#' @param ... Additional arguments to pass to ggplot2::scale_fill_manual.
#'
#' @return Returns a ggplot2 plot object for tissue images with optional
#' annotations.
#'
#' @details This function generates a ggplot2 plot for tissue images with
#' optional annotations, allowing customisation of various plot elements.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_tissueImgPlot
#'
#' @importFrom tidyterra geom_spatraster_rgb
#'
.int_tissueImgPlot <- function(image,
                               image_name,
                               limits_list,
                               sfe,
                               type,
                               annotate,
                               colours,
                               fill_args,
                               alpha,
                               ...) {
  ## Identify the RGB value range. It can be 0-1, 0-255 (8-bit), or 0-65536 (16-bit)
  ## There is a chance the 'range_max' slot to be filled with NaN. We need to
  ## account for that.
  max <- image@ptr$range_max[1]
  if (is.na(max)) {
    max <- max(image[,,1])
  }
  if (max <= 1) {
    max_col_value <- 1
  } else if (max <= 255) {
    max_col_value <- 255
  } else if (max <= 65536) {
    max_col_value <- 65536
  }

  ## To silence the 'coordinate system already present' warning as per:
  ## https://github.com/tidyverse/ggplot2/issues/2799
  cf <- coord_fixed()
  cf$default <- TRUE

  p <- ggplot() +
    tidyterra::geom_spatraster_rgb(data = image,
                                   max_col_value = max_col_value) +
    ggplot2::labs(subtitle = image_name) +
    ggplot2::lims(x = limits_list[[1]],
                  y = limits_list[[2]]) +
    cf +
    ggplot2::theme_void() +
    theme(plot.subtitle = ggplot2::element_text(hjust = 0.5))

  ## Put together the data to plot geometries. Choose spots or hexagons too.
  if (!type == "none") {
    ids <- image_name
    data <- .int_dataToPlotAnnot(sfe = sfe, type = type,
                                 ids = image_name, annotate = annotate)

    if (annotate) {
      ## Plot with annotation
      ## Fetch colours
      annot_cols <- .int_annotColours(sfe = sfe, colours = colours)
      ## Set some defaults if not provided
      fill_args <- .int_fillArgs(fill_args, annot_cols = annot_cols)
      ## Plot
      p <- p + ggplot2::geom_sf(data = data,
                                aes(geometry = geometry,
                                    fill = annotation),
                                alpha = alpha) +
        do.call(scale_fill_manual, c(list(...), fill_args)) +
        ggplot2::labs(fill = "Annotation")
    } else {
      ## Plot without annotation
      p <- p + ggplot2::geom_sf(data = data,
                                aes(geometry = geometry),
                                alpha = alpha)
    }
  }

  return(p)
}

# ---------------------------------------------------------------------------- #
# ------------------------------- DEPRECATED --------------------------------- #
# ---------------------------------------------------------------------------- #
#' Internal Function: .int_getImgLims2
#'
#' Extracts image limits for plotting based on spot locations for each sample.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#'
#' @return Returns a list of x and y limits for image plotting, organized by
#' sample IDs.
#'
#' @details This function extracts spot locations and sets x and y limits for
#' image plotting for each specified sample ID.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_getImgLims2
#'
.int_getImgLims2 <- function(sfe, sample_id = NULL) {
  ## This function is over-engineered. I forgot that the limits I get were the
  ## outer limits of the Visium capture area and not of the tissue's min max.
  ## I will keep it around in case I need to have different limits per sample.
  ## Fetch the required sample IDs
  sample_id <- .int_getImgSmplIDs(sfe = sfe, sample_id = sample_id)

  ## Find limits
  limits_list <- lapply(sample_id, function(i) {
    ## Extract the spot locations
    select_smpl <- colData(sfe)$sample_id == i
    spot_coords <- spatialCoords(sfe) %>% as.data.frame() %>% .[select_smpl,]

    ## Set limits
    xlim <- c(min(spot_coords$pxl_col_in_fullres) - 100,
              max(spot_coords$pxl_col_in_fullres) + 100)
    ylim <- c(min(spot_coords$pxl_row_in_fullres) - 100,
              max(spot_coords$pxl_row_in_fullres) + 100)

    return(c(xlim, ylim))
  })

  ## Add names
  names(limits_list) <- sample_id

  return(limits_list)
}
