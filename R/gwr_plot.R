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
plotGWR <- function(gwr, predictorVar = NULL) {
  ## Prepare the data for plotting
  ## Get the independent/predictor variable if not provided.
  ## Otherwise, a vector of variable names is expected!
  ## If NULL and there are multiple variables, only the first one is selected!
  if (is.null(predictorVar)) {
    mdl <- gwr$GW.arguments$formula
    predictorVar <- gsub(".+?~", "", mdl) %>%
      strsplit(., "[+]")
    predictorVar <- predictorVar[[1]][1]
  }

  predTVcol <- paste0(predictorVar, "_TV")
  gwr_sf <- gwr_toSF(gwr) %>%
    mutate(signif0 = factor(if_else(abs(.data[["Intercept_TV"]]) > 1.96,
                                    "Signif.",
                                    "Not-signif.")),
           signif1 = factor(if_else(abs(.data[[predTVcol]]) > 1.96,
                                    "Signif.",
                                    "Not-signif.")),
           b0_geom = if_else(signif0 == "Signif.",
                             geometry,
                             NA),
           b1_geom = if_else(signif1 == "Signif.",
                             geometry,
                             NA))

  ## Create table
  tab = round(summary(gwr$lm)$coefficients, 3)
  colnames(tab)[3:4] = c("t-value", "p-value")
  rownames(tab) = c("Intercept", "Covariate")
  tab = data.frame(tab)

  ## Plot
  p <- wrap_plots(.int_plotGWRCoeff(gwr_sf, 0, predictorVar),
                  .int_plotGWRCoeff(gwr_sf, 1, predictorVar),
                  ncol = 2)

  p <- p +
    annotation_custom(gridExtra::tableGrob(tab),
                      xmin = -3000, xmax = 200,
                      ymin = -1000, ymax = -300)

  p + plot_annotation(title = gwr$GW.arguments$formula,
                      theme = theme(plot.title = element_text(hjust = 0.5)))
}



#---------------------------------------------------------------------------- #
  #  ############# INTERNAL FUNCTIONS ASSOCIATED WITH MISC PLOTS ##############
# ---------------------------------------------------------------------------- #
#' Internal Function: Plot GWR Coefficients
#'
#' This internal function generates a spatial plot for the coefficients
#' obtained from a Geographically Weighted Regression (GWR) model. The function
#' highlights areas where the coefficients are statistically significant and
#' provides a visual representation of the intercept or a specified predictor
#' variable.
#'
#' @param gwr_sf A Simple Features (sf) object containing the spatial data and
#'   the GWR coefficients. This object should include columns for the
#'   coefficients and their significance, as well as the geometry information.
#' @param i An integer (0 or 1) indicating whether to plot the intercept (0) or
#'   the predictor variable (1) coefficients.
#' @inheritParams plotGWR
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
#' @rdname dot-int_plotGWRCoeff
#' @aliases .int_plotGWRCoeff
#'
#' @importFrom ggplot2 scale_fill_gradient2 alpha
#'
.int_plotGWRCoeff <- function(gwr_sf, i, predictorVar) {
  if (i == 0) {
    fill <-  "Intercept"
    b_geom <- "b0_geom"
    tit = expression(""*beta[0]*"")
  } else {
    fill  <- predictorVar
    b_geom <- "b1_geom"
    tit = expression(""*beta[1]*"")
  }

  ### Plot: Map Coefficients
  ggplot(gwr_sf) +
    geom_sf(aes(geometry = geometry,
                fill = get(fill)),
            colour = NA) +
    geom_sf(aes(geometry = get(b_geom)),
            fill = alpha("white", 0),
            colour = "black",
            linewidth = 0.15) +
    scale_fill_gradient2(high = "#880A1F",
                         mid = "#FAFBF7",
                         low = "#426E92",
                         midpoint = 0,
                         n.breaks = 7) +
    labs(x = "",
         y = "",
         fill = tit) +
    theme_void()
}
