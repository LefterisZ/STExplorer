#' Title: Visualise Global Spatial Autocorrelation result summary
#'
#' Description: This function plots a dot plot of global SA values for selected
#' genes.
#'
#' @param m_sfe An object of class 'SpatialFeatureExperiment' or
#'              'MultiSpatialFeatureExperiment' containing spatial data.
#' @param sample_id A character string specifying the sample ID for analysis.
#' @param genes A character vector containing a selection of genes. It must be
#'              either EnsgIDs OR gene names. A mix of both will throw an error.
#' @param statistic A character vector specifying which spatial autocorrelation
#'                  statistics to plot. Valid options are "moran", "getis", and
#'                  "geary". Default is c("moran", "getis", "geary").
#' @param ylab A character string indicating the label for the y-axis. It can
#'             be "gene_name" to use gene names or any other value to use
#'             Ensemble Gene IDs. Default is "gene_name".
#' @param moran_thresh A numeric vector of length 2 specifying the threshold
#'                     range for Moran's I values to be considered significant.
#'                     Default is c(-0.5, 0.5).
#' @param geary_thresh A numeric vector of length 2 specifying the threshold
#'                     range for Geary's C values to be considered significant.
#'                     Default is c(0.5, 1.5).
#' @param getis_thresh A numeric value specifying the threshold for Getis-Ord's
#'                     G statistic to be considered significant. Default is 2.
#' @param ... Additional arguments passed on to `geom_point()` for customising
#'            the point aesthetics in the ggplot.
#'
#' @return A ggplot object representing the dot plot of global spatial
#'         autocorrelation statistics for the specified genes across the given
#'         statistics.
#'
#' @details
#' The function creates a faceted dot plot where each facet corresponds to one
#' of the specified spatial autocorrelation statistics. Points represent genes,
#' and their position and colour indicate the statistic value. Point size shows
#' the significance of the statistic, while the shape indicates whether the
#' statistic suggests positive or negative spatial autocorrelation.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial autocorrelation visualisation global
#'
#' @rdname plotSA_globalSum
#'
#' @importFrom ggplot2 theme_bw
#'
#' @examples
#' plotSA_globalSum(m_sfe, "sample1", c("Gene1", "Gene2"), c("moran", "getis"))
#'
#' @export
plotSA_globalSum <- function(m_sfe,
                             sample_id = NULL,
                             genes,
                             statistic = c("moran", "getis", "geary"),
                             ylab = "gene_name",
                             moran_thresh = c(-0.5, 0.5),
                             geary_thresh = c(0.5, 1.5),
                             getis_thresh = 2,
                             ...) {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Validate genes input
  gs_to_plot <- .int_matchNamesToEnsgID(sfe, genes)

  ## Select columns to plot
  stat_to_plot <- .int_getSAStatColumns(colnames(rowData(sfe)),
                                        statistic = statistic,
                                        ylab = ylab)

  ## Prepare data to plot
  data <- .int_prepareGlobalSAPlotData(rowData(sfe)[gs_to_plot, stat_to_plot],
                                       moran_thresh = moran_thresh,
                                       geary_thresh = geary_thresh,
                                       getis_thresh = getis_thresh)

  ## Prepare plot elements
  if (ylab == "gene_name") {
    y_lab <- "Gene Names"
  } else {
    y_lab <- "Ensemble Gene IDs"
  }

  ## Plot
  ggplot(data = data,
         aes(x = stat, y = gene_name,
             colour = stat,
             shape = PosNegSA,
             size = Pval)) +
    geom_point(...) +
    facet_wrap(~stat_name, scales = "free_x") +
    labs(y = y_lab,
         x = "SA Statistic",
         size = "p-value",
         colour = "SA\nStatistic",
         shape = "+/- SA") +
    theme_bw()
}


#' Title: Visualise Local Spatial Autocorrelation
#'
#' Description: Generates a map visualizing local spatial autocorrelation
#' values for a given feature in the input spatial feature experiment.
#'
#' @param m_sfe An object of class 'SpatialFeatureExperiment' or
#'              'MultiSpatialFeatureExperiment' containing spatial data.
#' @param sample_id A character string specifying the sample ID for analysis.
#' @param feature A character string indicating the feature of interest.
#'                It can be either ENSGene IDs or Gene Names.
#' @param statistic A character string specifying the type of spatial
#'                  autocorrelation results ("moran", "geary", "getis") to look
#'                  for in the m_sfe object for the specified sample.
#' @param test A character string specifying the statistical test method
#'              ("z-score" or "permutation") to look for in the m_sfe object
#'              for the specified sample.
#' @param pVal A numeric value representing the significance threshold for
#'             plotting. Ignored for Geary's C, where permutation is suggested.
#'             See details for more information.
#' @param type A character string specifying the type of geometries ("spot" or
#'             "hex") to use in the plot.
#' @param viridis_col A character string specifying the colour palette for
#'                    the autocorrelation values. Currently supporting only the
#'                    viridis palettes.
#' @param signif_col A character string specifying the colour for the outline
#'                   of significant locations in the plot. Used only when the
#'                   `locations` argument is set to `"all"`. It helps to
#'                   identify which SA values are also statistically
#'                   significant.
#' @param title A character string ("ENSGID" or "name") specifying whether the
#'              title should show the ENSGene ID or the gene name. Default is
#'              "name".
#' @param locations A character string ("all" or "significant") specifying
#'                  which locations to include in the plot. If `"all"` is
#'                  selected, then all spots are included and the locations
#'                  where the local SA is significant has a differently
#'                  coloured outline as set by the `signif_col` argument. If
#'                  `"significant` is selected, then only the locations with
#'                  statistically significant SA values are plotted while the
#'                  rest of the locations are greyed out.
#'
#' @return A ggplot object visualizing local spatial autocorrelation.
#'
#' @details Creates an informative plot of local spatial autocorrelation values.
#'
#' Internally this function checks the compatibility between selected statistic
#' and statistical significance method. Geary's C in particular returns no
#' Z-score statistic. As a result, when someone is trying to plot Geary's C the
#' `pVal` argument is not used. A message will be printed to
#' inform for this and also to suggest using `permutation` in the `test`
#' argument after running the `gearyLocalCPerm` function first.
#'
#' @seealso See also the `plotSA_localClust` function for clustering analysis.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial, autocorrelation, visualisation
#'
#' @rdname plotSA_local
#'
#' @importFrom ggplot2 ggplot scale_fill_viridis_c geom_sf theme_void theme labs
#'
#' @examples
#' plotSA_local(m_sfe, "sample1", "gene_expression", "moran", "z-score",
#' 0.05, "spot")
#'
#' @export
plotSA_local <- function(m_sfe,
                         sample_id = NULL,
                         feature,
                         statistic = c("moran", "geary", "getis"),
                         test = c("z-score", "permutation"),
                         pVal = 0.05,
                         type = c("spot", "hex"),
                         viridis_col = "inferno",
                         signif_col = "#E0E0E0",
                         title = c("ENSGID", "name"),
                         locations = "all") {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check valid input arguments
  statistic <- match.arg(statistic)
  test <- match.arg(test)
  type <- match.arg(type)
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
  }
  if (missing(title)) {
    title <- "name"
  }

  ## Validate genes input
  feature <- .int_matchNamesToEnsgID(sfe, feature)

  ## Construct plot's title
  title <- .int_getTitle(sfe = sfe, feature = feature, title = title)

  ## Get the `localResults` name for the specific statistic
  localRes <- .int_getLocalResName(statistic = statistic, test = test)

  ## Check if `localGearyC` is being used. If yes, set `pVal` to NULL
  pVal <- .int_checkSACompatible(localRes = localRes, pVal = pVal)

  ## Prepare the data to plot
  data <- .int_prepareLocalSAPlotData(sfe = sfe,
                                      feature = feature,
                                      localRes = localRes,
                                      pVal = pVal,
                                      type = type,
                                      plot = "localStat")

  ## Construct plot's fill label
  fill <- .int_getFillLabel(statistic = statistic, plot = "localStat")

  ## Workaround to suppress the 'Coordinate system already present' message
  ## from ggplot2 based on this GitHub issue:
  ## https://github.com/tidyverse/ggplot2/issues/2799
  cf <- coord_fixed()
  cf$default <- TRUE

  if (statistic == "geary") {
    direction <- -1
  } else {
    direction <- 1
  }

  ## Create the ggplot object
  p <- ggplot2::ggplot(data = data) +
    ggplot2::scale_fill_viridis_c(option = viridis_col,
                                  direction = direction) +
    cf +
    #ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(fill = fill,
                  title = title)

  ## Add geom_sf
  ## - for an unknown reason, when the geom_sf() is added inside the above it
  ##   throws an error. Temporary solution until fixed:
  if (locations == "all") {
    p <- p + geom_sf(aes(geometry = geometry, fill = Stat),
                     colour = NA)
    p <- p + geom_sf(aes(geometry = geom_b),
                     fill = NA, colour = signif_col, linewidth = 0.25)
  } else if (locations == "significant") {
    p <- p + geom_sf(aes(geometry = geometry),
                     colour = NA, fill = "#E0E0E0", alpha = 0.5)
    p <- p + geom_sf(aes(geometry = geom_b, fill = Stat),
                     colour = NA)
  }

  ## Return plot
  return(p)
}


#' #' Title: Visualise Clusters in Local Spatial Autocorrelation
#'
#' Description: Generates a map visualising clusters in local spatial
#' autocorrelation values for a given feature in the input spatial feature
#' experiment.
#'
#' @param m_sfe An object of class 'SpatialFeatureExperiment' or
#'              'MultiSpatialFeatureExperiment' containing spatial data.
#' @param sample_id A character string specifying the sample ID for analysis.
#' @param feature A character string indicating the feature of interest.
#'                Currently supports only ENSGene IDs.
#' @param statistic A character string specifying the type of spatial
#'                  autocorrelation results ("moran", "geary", "getis") to look
#'                  for in the m_sfe object for the specified sample.
#' @param test A character string specifying the statistical test method
#'              ("z-score" or "permutation") to look for in the m_sfe object
#'              for the specified sample.
#' @param pVal A numeric value representing the significance threshold for
#'             plotting. Ignored for Geary's C, where permutation is suggested.
#' @param type A character string specifying the type of geometries ("spot" or
#'             "hex") to use in the plot.
#' @param title A character string ("ENSGID" or "name") specifying whether the
#'              title should show the ENSGene ID or the gene name.
#' @param clust_col A named character vector of colours to be used for
#'                  plotting. The names must match the cluster names. Look at
#'                  details below for more information about what these names
#'                  should be. Default is `NULL` and the default colours are
#'                  applied.
#' @param signif_col A character string specifying the colour for the outline
#'                   of significant locations in the plot. Used only when the
#'                   `locations` argument is set to `"all"`. It helps to
#'                   identify which SA values are also statistically
#'                   significant.
#' @param locations A character string ("all" or "significant") specifying
#'                  which locations to include in the plot. If `"all"` is
#'                  selected, then all spots are included and the locations
#'                  where the local SA is significant has a differently
#'                  coloured outline as set by the `signif_col` argument. If
#'                  `"significant` is selected, then only the locations with
#'                  statistically significant SA values are plotted while the
#'                  rest of the locations are greyed out.
#'
#' @return A ggplot object visualising clusters in local spatial
#'         autocorrelation.
#'
#' @details Creates an informative plot of local spatial autocorrelation
#'          values.
#'          The named vector with colours provided in the `clust_col` argument
#'          should include the below names:
#'          "Low", "High", "low-Low", "High-Low", "Low-High", "High-High",
#'          "Not Signif."
#'          If you don't provide a vector with the exact names then this
#'          function will produce an error.
#'
#' @seealso See also the 'plotSA_local' function for local autocorrelation
#'          visualization.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords spatial, autocorrelation, clustering, visualisation
#'
#' @rdname plotSA_localClust
#'
#' @importFrom ggplot2 ggplot scale_fill_viridis_c geom_sf theme_void theme labs
#'
#' @examples
#' plotSA_localClust(m_sfe, "sample1", "gene_expression", "moran", "z-score",
#' 0.05, "spot")
#'
#' @export
plotSA_localClust <- function(m_sfe,
                              sample_id = NULL,
                              feature,
                              statistic = c("moran", "geary", "getis"),
                              test = c("z-score", "permutation"),
                              pVal = 0.05,
                              type = c("spot", "hex"),
                              title = c("ENSGID", "name"),
                              clust_col = NULL,
                              signif_col = "#E0E0E0",
                              locations = "all") {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check valid input arguments
  statistic <- match.arg(statistic)
  test <- match.arg(test)
  type <- match.arg(type)
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
  }
  if (!is.null(clust_col)) {
    if (is.null(names(clust_col))) {
      stop("The colour vector provided to the 'clust_col' argument MUST be ",
           "a NAMED vector.\nFor the names have a look at the help page:\n",
           "?plotSA_localClust")
    }
  }
  if (missing(title)) {
    title <- "name"
  }

  ## Validate genes input
  feature <- .int_matchNamesToEnsgID(sfe, feature)

  ## Construct plot's title
  title <- .int_getTitle(sfe = sfe, feature = feature, title = title)

  ## Get the `localResults` name for the specific statistic
  localRes <- .int_getLocalResName(statistic = statistic, test = test)

  ## Check if `localGearyC` is being used. If yes, set `pVal` to NULL
  pVal <- .int_checkSACompatible(localRes = localRes, pVal = pVal)

  ## Prepare the data to plot
  data <- .int_prepareLocalSAPlotData(sfe = sfe,
                                      feature = feature,
                                      localRes = localRes,
                                      pVal = pVal,
                                      type = type,
                                      plot = "cluster")

  ## Construct plot's fill label
  fill <- .int_getFillLabel(statistic = statistic, plot = "cluster")

  ## Get colours
  colours <- .int_getSAClustColours(clust_col = clust_col,
                                    clusters = levels(data$Clust))

  ## Workaround to suppress the 'Coordinate system already present' message
  ## from ggplot2 based on this GitHub issue:
  ## https://github.com/tidyverse/ggplot2/issues/2799
  cf <- coord_fixed()
  cf$default <- TRUE

  ## Create the ggplot object
  p <- ggplot2::ggplot(data = data) +
    ggplot2::scale_fill_manual(values = colours) +
    cf +
    #ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(fill = fill,
                  title = title)

  ## Add geom_sf
  ## - for an unknown reason, when the geom_sf() is added inside the above it
  ##   throws an error. Temporary solution until fixed:
  if (locations == "all") {
    p <- p + geom_sf(aes(geometry = geometry, fill = Clust),
                     colour = NA)
    p <- p + geom_sf(aes(geometry = geom_b),
                     fill = NA, colour = signif_col, linewidth = 0.25)
  } else if (locations == "significant") {
    p <- p + geom_sf(aes(geometry = geometry),
                     colour = NA, fill = "#E0E0E0", alpha = 0.5)
    p <- p + geom_sf(aes(geometry = geom_b, fill = Clust),
                     colour = NA)
  }

  ## Return plot
  return(p)
}


# ---------------------------------------------------------------------------- #
#  ######### INTERNAL FUNCTIONS ASSOCIATED WITH SA PLOTS  (C, G, I) #########
# ---------------------------------------------------------------------------- #
#' Internal Function: Select global SA stat columns to plot
#'
#' This function is fetching the required columns from the SFE's rowData.
#'
#' @param colnames The rowData column names from an SFE object.
#' @inheritParams plotSA_globalSum
#'
#' @return A character string representing the name of the local results.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getSAStatColumns
#'
.int_getSAStatColumns <- function(colnames, statistic, ylab) {
  if ("moran" %in% statistic) {
    col_moran <- grep("moran", colnames)
  } else {
    col_moran <- NULL
  }
  if ("geary" %in% statistic) {
    col_geary <- grep("geary", colnames)
  } else {
    col_geary <- NULL
  }
  if ("getis" %in% statistic) {
    col_getis <- grep("getis", colnames)
  } else {
    col_getis <- NULL
  }

  if (ylab == "gene_name") {
    col_g_name <- grep("gene_name", colnames)
  } else {
    col_g_name <- NULL
  }

  stat_to_plot <- c(col_g_name, col_moran, col_geary, col_getis)

  return(stat_to_plot)
}


#' Internal Function: Prepare Global SA summary data from plotting
#'
#' @param data A subset of SFE's rowData
#' @inheritParams plotSA_globalSum
#'
#' @return A long format data frame.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_prepareGlobalSAPlotData
#'
#' @aliases .int_prepareGlobalSAPlotData
#'
.int_prepareGlobalSAPlotData <- function(data,
                                         moran_thresh,
                                         geary_thresh,
                                         getis_thresh) {
  data <- data %>%
    as.data.frame() %>%
    rownames_to_column(var = "g_id")

  # Remove from colnames "Perm" and "Test"
  colnames(data) <- gsub("Perm|Test", "", colnames(data))

  data <- data %>%
    pivot_longer(cols = -any_of(c("g_id","gene_name")),
                 names_to = c("stat_name", ".value"),
                 names_sep = "_") %>%
    dplyr::mutate(stat_name = dplyr::case_match(stat_name,
                                                "moranI" ~ "Moran's I",
                                                "gearyC" ~ "Geary's C",
                                                "getisG" ~ "Getis & Ord's G"),
                  PosNegSA = dplyr::case_when(stat_name == "Moran's I" &
                                                stat <= moran_thresh[1] ~ "-",
                                              stat_name == "Moran's I" &
                                                stat >= moran_thresh[2] ~ "+",
                                              stat_name == "Geary's C" &
                                                stat <= geary_thresh[1] ~ "+",
                                              stat_name == "Geary's C" &
                                                stat >= geary_thresh[2] ~ "-",
                                              stat_name == "Getis & Ord's G" &
                                                stat >= getis_thresh ~ "+",
                                              .default = "No"))
}


#' Internal Function: Get localResults Name for Spatial Autocorrelation
#'
#' @param statistic A character string specifying the type of spatial
#'                  autocorrelation used ("moran", "geary", "getis").
#'
#' @param test A character string specifying the statistical test method used
#'             ("z-score" or "permutation").
#'
#' @return A character string representing the name of the local results.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getLocalResName
#'
#' @aliases .int_getLocalResName
#'
.int_getLocalResName <- function(statistic, test) {
  ## Set test part of the name
  if (test == "z-score") {
    test <- NULL
  } else {
    test <- "Perm"
  }

  ## Set statistic's name part of the name
  if (statistic == "moran") {
    stat <- "MoranI"
  } else if (statistic == "geary") {
    stat <- "GearyC"
  } else if (statistic == "getis") {
    stat <- "GetisG"
  }

  paste0("local",stat,test)
}


#' Internal Function: Check compatibility between selected statistic and
#' statistical significance method
#'
#' This internal function checks the compatibility between selected statistic
#' and statistical significance method. Geary's C in particular returns no
#' Z-score statistic. As a result, when someone is plotting Geary's C the
#' `pVal` argument is not used. A message will be printed to
#' inform for this and also to suggest using `permutation` instead.
#'
#' @param localRes A character string with the `localResult` name.
#' @param pVal Numeric. The p-value.
#'
#' @details
#' This function is intended for internal use to validate the compatibility
#' of the input parameters for plotting local spatial autocorrelation
#' statistics. It checks if `geary` is provided with `z-score`, if yes, then
#' sets the `pVal` argument to `NULL`. Otherwise returns the `pVal` argument as
#' is.
#'
#' @return This function returns a value for the `pVal` argument.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_checkSACompatible
#'
#' @aliases .int_checkSACompatible
.int_checkSACompatible <- function(localRes, pVal) {
  if (localRes == "localGearyC") {
    warning("The local Geary's C does not include a 'z-score' calculation to ",
            "infer significance of the results. \n",
            "1.  As a result, the pVal argument is not used. \n",
            "2.  If you want to filter for the 'interesting' locations, we ",
            "suggest using the 'gearyLocalCPerm' function to calculate \n",
            "    statistical significance using permutations. \n",
            "3.  Then, privide 'permutation' in the 'signif' argument.")
    pVal <- NULL
  }

  return(pVal)
}


#' Internal Function: Prepare Data for Spatial Autocorrelation Plot
#'
#' @param sfe An object of class 'SpatialFeatureExperiment' containing spatial
#'            data.
#' @param feature A character string indicating the feature of interest.
#' @param localRes A character string with the name of the local results.
#' @param pVal A numeric value representing the statistical significance
#'             threshold.
#' @param type A character string specifying the type of plot ("spot" or "hex").
#' @param plot A character string specifying the type of plot ("localStat" or
#'             "cluster").
#'
#' @importFrom dplyr case_when
#' @importFrom SpatialFeatureExperiment localResult
#'
#' @return A data frame containing the data for spatial autocorrelation plot.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_prepareLocalSAPlotData
#'
#' @aliases .int_prepareLocalSAPlotData
#'
.int_prepareLocalSAPlotData <- function(sfe,
                                   feature,
                                   localRes,
                                   pVal,
                                   type,
                                   plot = c("localStat","cluster")) {
  ## Select appropriate geometries
  if (type == "spot") {
    geoms <- "spotPoly"
  } else if (type == "hex") {
    geoms <- "spotHex"
  }

  ## Export local SA results
  dt <- data.frame(Stat = localResult(x = sfe,
                                      type = localRes,
                                      feature = feature)[[1]],
                   p.value = localResult(x = sfe,
                                         type = localRes,
                                         feature = feature)[[2]],
                   geometry = colGeometries(sfe)[[geoms]][["geometry"]]) %>%
          mutate(signif = as.factor(case_when(p.value < pVal ~ "Signif.",
                                              p.value >= pVal ~ "Not Signif.")),
                 geom_b = if_else(signif == "Signif.", .data$geometry, NA))

  if (plot == "cluster") {
    dt <- dt %>%
      mutate(Clust = localResult(x = sfe,
                           type = localRes,
                           feature = feature)[[3]])

  }
  return(dt)
}


#' Internal Function: Get Fill Label for Spatial Autocorrelation Plot
#'
#' @param statistic A character string specifying the type of spatial
#'                  autocorrelation used ("moran", "geary", "getis").
#'
#' @param plot A character string specifying the type of plot ("localStat" or
#'             "cluster").
#'
#' @return A character string representing the fill label for the plot.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getFillLabel
#'
#' @aliases .int_getFillLabel
#'
.int_getFillLabel <- function(statistic,
                              plot = c("localStat", "cluster")) {
  if (statistic == "moran") {
    fill <- "Local\nMoran's\nIi"
  } else if (statistic == "geary") {
    fill <- "Local\nGeary's\nCi"
  } else if (statistic == "getis") {
    fill <- "Local\nGetis's\nGi"
  }

  ## If cluster add "\nClusters" at the end
  if (plot == "cluster") {
    fill <- paste0(fill, "\nClusters")
  }

  return(fill)
}


#' Internal Function: Get Title for Spatial Autocorrelation Plot
#'
#' @param sfe An object of class 'SpatialFeatureExperiment'
#'            containing spatial data.
#'
#' @param feature A character string indicating the feature of interest.
#'
#' @param title A character string specifying the title ("ENSGID" or
#'              "name").
#'
#' @return A character string representing the title for the plot.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getTitle
#'
#' @aliases .int_getTitle
#'
.int_getTitle <- function(sfe,
                          feature,
                          title) {
  if (title == "ENSGID") {
    out <- feature
  } else if (title == "name") {
    out <- rowData(sfe)[rowData(sfe)$id == feature, "gene_name"]
  }

  return(out)
}


#' Internal Function: Get Colours for Spatial Autocorrelation Clusters
#'
#' @param clust_col A named character vector of colours for plotting clusters.
#'                  The names must match cluster names.
#'
#' @param clusters A character vector of cluster names.
#'
#' @return A character vector of colours for spatial autocorrelation clusters.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_getSAClustColours
#'
#' @aliases .int_getSAClustColours
#'
.int_getSAClustColours <- function(clust_col, clusters) {
  cols_list <- list("Low" = "#0066ff",
                    "High" = "#ff0000",
                    "Low-Low" = "#0066ff",
                    "High-Low" = "#ffb3b3",
                    "Low-High" = "#99c2ff",
                    "High-High" = "#ff0000",
                    "Not Signif." = "#E0E0E0",
                    "Other Positive" = "#99c2ff",
                    "Negative" = "#ffb3b3")

  if (!is.null(clust_col)) {
    cols_list[names(clust_col)] <- clust_col
  }

  out <- cols_list[clusters]
  out <- unlist(out)

  return(out)

}
