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

  gwr_sf <- st_as_sf(gwr$SDF) %>%
    mutate(signif0 = factor(if_else(abs(.data[["Intercept_TV"]]) > 1.96,
                                    "Signif.",
                                    "Not-signif.")),
           signif1 = factor(if_else(abs(.data[[paste0(predictorVar, "_TV")]]) > 1.96,
                                    "Signif.",
                                    "Not-signif.")),
           b0_geom = if_else(signif0 == "Signif.",
                             geometry,
                             NA),
           b1_geom = if_else(signif1 == "Signif.",
                             geometry,
                             NA))

  for (i in c(0,1)) {
    tab = round(summary(gwr$lm)$coefficients, 3)
    colnames(tab)[3:4] = c("t-value", "p-value")
    rownames(tab) = c("Intercept", "Covariate")

    tab = data.frame(tab)
    tab1 = tab[1,]
    tab2 = tab[2,]

    if (i == 0) {
      fill <-  "Intercept"
      b_geom <- "b0_geom"
      tit = expression(""*beta[0]*"")
      table = tab1
    } else {
      fill  <- predictorVar
      b_geom <- "b1_geom"
      tit = expression(""*beta[1]*"")
      table = tab2
    }

    ### Plot: Map Coefficients that flip ----
    ggplot(gwr_sf) +
       geom_sf(aes(geometry = geometry,
                  fill = get(fill)),
               colour = NA) +
       geom_sf(aes(geometry = get(b_geom)),
               fill = alpha("white", 0),
               colour = "black",
               linewidth = 0.5) +
       scale_fill_gradient2(high = "#880A1F",
                            mid = "#FAFBF7",
                            low = "#426E92",
                            midpoint = 0,
                            n.breaks = 7) +
       labs(x = "",
            y = "",
            fill = tit) +
       annotation_custom(tableGrob(table),
                         xmin = 150, xmax = 450,
                         ymin = -570, ymax = -500) +
       theme_void()
  }
}

