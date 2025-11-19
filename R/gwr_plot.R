#' Plot Geographically Weighted Regression (GWR) Results
#'
#' This function generates a plot that visualises the results of a
#' Geographically Weighted Regression (GWR) model. The function plots
#' the coefficients of the GWR model, highlighting significant areas
#' based on t-values, and includes a summary table of the model's coefficients.
#'
#' @param gwr A GWR object returned by the GWR fitting function `gwrSTE`.
#' @param predictorVar A character string specifying the name of the predictor
#'   variable to plot. If `NULL`, the function will automatically select the
#'   first predictor variable from the model formula. If the model contains
#'   multiple variables as predictors, only the first one will be used.
#' @param globalModel TRUE or FALSE. Indicates whether a table with the global
#'   model's summary statistics has to be printed.
#' @param title A character string for the plot's title. Defaults to NULL where
#'   a title made from the formula is used.
#'
#' @details
#' The generated plot shows the spatial distribution of the GWR coefficients
#' and their significance. Areas where the coefficients are significant are
#' highlighted. The plot also includes a summary table of the global model's
#' intercept and covariate coefficients, t-values, and p-values to allow the
#' user to compare results with GWR (the mapped coefficients).
#'
#' @return A `ggplot` object containing the plot of GWR coefficients and a
#'   summary table.
#'
#' @importFrom ggplot2 annotation_custom
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom grid unit
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname plotGWR
#'
#' @examples
#' \dontrun{
#' pltGWR(gwr_object)
#' }
#'
#'
#' @export
plotGWR <- function(gwr,
                    predictorVar = NULL,
                    globalModel = FALSE,
                    title = NULL) {
  ## Prepare the data for plotting
  ## Get the independent/predictor variable if not provided.
  ## Otherwise, a vector of variable names is expected!
  if (is.null(predictorVar)) {
    mdl <- gwr$GW.arguments$formula
    predictorVar <- gsub(".+?~", "", mdl) %>%
      strsplit(., "[+]")
    predictorVar <- predictorVar[[1]]
  }

  ## Get all variables to plot
  vars_to_plot <- c("Intercept", predictorVar)

  ## Create table
  if (globalModel) {
    tab = round(summary(gwr$lm)$coefficients, 3)
    colnames(tab)[3:4] = c("t-value", "p-value")
    if (nrow(tab) == 2) {
      rownames(tab) = c("Intercept", "Covariate")
    } else {
      rownames(tab) = c("Intercept", paste0("Covariate", 1:(nrow(tab) - 1)))
    }
    tab = data.frame(tab)
  }

  ## Prepare title
  if (is.null(title)) {
    title <- gwr$GW.arguments$formula
  } else if (!is.character(title)) {
    stop("The title argument must be either left to NULL or provided with a ",
         "character string.")
  }

  ## Plot
  plot_list <- lapply(seq_along(vars_to_plot),
                      function(x) {
                        Var = vars_to_plot[x]
                        dt <- STExplorer:::.int_getGWRdataToPlot(gwr = gwr, Var = Var)
                        p <- .int_plotGWRCoeff(gwr_sf = dt, i = x, Var = Var)
                        return(p)
                      })

  p <- wrap_plots(plot_list, ncol = 2)

  if (globalModel) {
    p <- p |
      ggplot() +
      annotation_custom(gridExtra::tableGrob(tab),
                        xmin = -Inf, xmax = Inf,
                        ymin = -Inf, ymax = Inf) +
      theme_void() +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"))
  }

  p + plot_annotation(title = title,
                      theme = theme(plot.title = element_text(hjust = 0.5)))
}



#' Plot Local R-squared from GWR Model
#'
#' This function generates a plot of the local R-squared values from a
#' Geographically Weighted Regression (GWR) model. The plot provides a
#' spatial visualization of the goodness of fit across the study area.
#'
#' @param gwr A GWR object returned by the GWR fitting function `gwrSTE`.
#'
#' @return A `ggplot` object representing the spatial distribution of the
#'   local R-squared values.
#'
#' @importFrom ggplot2 geom_sf scale_fill_viridis_c labs theme_void
#' @importFrom dplyr select
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname plotGWR_R2
#'
#' @examples
#' \dontrun{
#' plotGWR_R2(gwr_object)
#' }
#'
#' @export
plotGWR_R2 <- function(gwr) {
  ## Extract R-squared data
  gwr_sf <- gwr_toSF(gwr) %>%
    select(all_of(c("Local_R2", "geometry")))

  ## Plot R-squared
  ggplot(gwr_sf) +
    geom_sf(aes(geometry = geometry, fill = Local_R2),
            colour = NA) +
    scale_fill_viridis_c() +
    labs(fill = expression("R"^2)) +
    theme_void()
}

#---------------------------------------------------------------------------- #
  #  ############# INTERNAL FUNCTIONS ASSOCIATED WITH MISC PLOTS ##############
# ---------------------------------------------------------------------------- #
#' Internal Function: Prepare GWR Data for Plotting
#'
#' This internal function prepares the data from a Geographically Weighted
#' Regression (GWR) model for plotting. It extracts and formats the coefficients
#' and their significance for visualization.
#'
#' @param gwr A GWR object containing the fitted GWR model results.
#' @param Var A character string specifying the variable whose coefficients
#'   are to be plotted.
#'
#' @return A Simple Features (sf) object containing the coefficients, their
#'   significance, and the geometry for plotting.
#'
#' @details
#' The function adds a column `signifCoef` to indicate the significance of
#' the coefficients based on the t-value threshold of 1.96. It also prepares
#' a `geom_signif` column to facilitate the plotting of significant areas.
#'
#' This function is used internally by the `plotGWR` function to prepare
#' spatial data for plotting.
#'
#' @seealso
#' \code{\link{plotGWR}}, \code{\link{gwr_toSF}}
#'
#' @keywords internal
#'
#' @importFrom dplyr mutate select if_else
#' @importFrom rlang sym
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_getGWRdataToPlot
#' @aliases .int_getGWRdataToPlot
#'
.int_getGWRdataToPlot <- function(gwr, Var) {
  ## Fetch predictor TV column
  VarTVcol <- paste0(Var, "_TV")

  ## Prepare the data
  gwr_sf <- gwr_toSF(gwr) %>%
    mutate(signifCoef = factor(if_else(abs(.data[[VarTVcol]]) > 1.96,
                                       "Signif.",
                                       "Not-signif.")),
           geom_signif = if_else(signifCoef == "Signif.",
                               geometry,
                               NA)) %>%
    select(all_of(c(Var, "geometry", "geom_signif")))

  return(gwr_sf)
}

#' Internal Function: Plot GWR Coefficients (Visium HD-safe)
#'
#' This internal function generates a spatial plot for the coefficients
#' obtained from a Geographically Weighted Regression (GWR) model. The function
#' highlights areas where the coefficients are statistically significant and
#' provides a visual representation of the intercept or a specified predictor
#' variable.
#'
#' @param gwr_sf sf with columns: Var, geometry, geom_signif (as created by .int_getGWRdataToPlot)
#' @param i integer (1 = intercept)
#' @param Var character, name of coefficient column to plot
#' @param use_alpha logical, if TRUE use alpha to show significance
#' @param alpha_range numeric length-2, alpha range for Not-sig → Sig
#' @param show_stipple logical, add light point stipple for significant loci
#' @param stipple_size numeric, point size for stipple (small for Visium HD)
#' @param stipple_alpha numeric in [0,1], opacity for stipple
#'
#' @return A `ggplot` object representing the spatial distribution of the
#'   specified GWR coefficient, with significant areas highlighted.
#'
#' @details
#' This function is an internal utility used to plot the spatial distribution
#' of GWR model coefficients. Depending on the value of the `i` parameter, the
#' function plots either the intercept or the specified predictor variable. The
#' plot uses a colour gradient to represent the magnitude of the coefficients
#' and overlays the significant areas with a thin black outline.
#'
#' The function is designed to be used within the `plotGWR` function and is not
#' intended for direct use by end-users.
#'
#' @seealso
#' \code{\link{plotGWR}}, \code{\link{gwr_toSF}}
#'
#' @keywords internal GWR coefficients spatial plot ggplot
#'
#' @author Eleftherios (Lefteris) Zormpas
#' @rdname dot-int_plotGWRCoeff
#' @aliases .int_plotGWRCoeff
#'
#' @importFrom ggplot2 scale_fill_gradient2 alpha
#'
.int_plotGWRCoeff <- function(gwr_sf, i, Var,
                              use_alpha     = TRUE,
                              alpha_range   = c(0.35, 1),
                              show_stipple  = FALSE,
                              stipple_size  = 0.05,
                              stipple_alpha = 0.25) {

  # Title (β0 for intercept; β_j otherwise)
  if (i == 1) {
    tit <- expression(""*beta[0]*"")
  } else {
    j <- i - 1
    tit <- bquote(beta[.(j)])
  }

  # Derive a lightweight significance flag from 'geom_signif'
  # ('.int_getGWRdataToPlot' already put the significant geometry there, NA otherwise)
  if ("geom_signif" %in% names(gwr_sf)) {
    # sfc supports is.na(); TRUE where significant geometry exists
    sig_flag <- !is.na(gwr_sf$geom_signif)
  } else {
    sig_flag <- rep(TRUE, nrow(gwr_sf))  # fallback: no gating
  }
  gwr_sf$..alpha_sig <- if (use_alpha) ifelse(sig_flag, 1, 0) else 1

  # Base coefficient layer, NO outlines
  p <- ggplot(gwr_sf) +
    geom_sf(aes(geometry = geometry,
                fill     = !!rlang::sym(Var),
                alpha    = ..alpha_sig),
            colour = NA) +
    scale_alpha_continuous(range = alpha_range, guide = "none") +
    scale_fill_gradient2(high     = "#880A1F",
                         mid      = "#FAFBF7",
                         low      = "#426E92",
                         midpoint = 0,
                         n.breaks = 7) +
    labs(x = NULL, y = NULL, fill = tit) +
    theme_void()

  # Optional: add a very light stipple (points) at significant centroids
  if (show_stipple && any(sig_flag, na.rm = TRUE) && "geom_signif" %in% names(gwr_sf)) {
    # Convert significant polygons to points (centroids/surface points)
    # Keep only rows that are flagged significant
    sig_rows <- which(sig_flag %in% TRUE)
    if (length(sig_rows)) {
      # st_point_on_surface avoids invalid centroids for concave polygons
      sig_pts <- try(sf::st_point_on_surface(gwr_sf$geometry[sig_rows]), silent = TRUE)
      if (!inherits(sig_pts, "try-error")) {
        pts_sf <- sf::st_sf(sig_pt = sig_pts)
        p <- p + geom_sf(data = pts_sf,
                         aes(geometry = sig_pt),
                         inherit.aes = FALSE,
                         colour = "black",
                         size   = stipple_size,
                         alpha  = stipple_alpha)
      }
    }
  }

  p
}
