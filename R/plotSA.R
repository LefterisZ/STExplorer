#' Map Local Spatial Autocorrelation values
#'
#' This function generates a map of the local spatial autocorrelation values.
#'
#'
#' @export
plotSA_local <- function(m_sfe,
                         sample_id,
                         feature,
                         statistic = c("moran", "geary", "getis"),
                         test = c("z-score", "permutation"),
                         pVal = 0.05,
                         type = c("spot", "hex"),
                         viridis_col = "inferno",
                         signif_col = "#E0E0E0") {
  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check valid input arguments
  statistic <- match.arg(statistic)
  test <- match.arg(test)
  type <- match.arg(type)
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
  }

  ## Get feature suffixes
  feature <- .int_getSAFeatSuffix(feature = feature, statistic = statistic)

  ## Get the `localResults` name for the specific statistic
  localRes <- .int_getLocalResName(statistic = statistic, test = test)

  ## Check if `localGearyC` is being used. If yes, set `pVal` to NULL
  pVal <- .int_checkSACompatible(localRes = localRes, pVal = pVal)

  ## Prepare the data to plot
  data <- .int_prepareSAPlotData(sfe = sfe,
                                 feature = feature,
                                 localRes = localRes,
                                 pVal = pVal,
                                 type = type,
                                 plot = "localStat")

  ## Construct plot's fill label
  fill <- .int_getFillLabel(statistic = statistic, plot = "localStat")

  ## Create the ggplot object
  p <- ggplot2::ggplot(data = data) +
    ggplot2::scale_fill_viridis_c(option = viridis_col) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(fill = fill)

  ## Add geom_sf
  ## - for an unknown reason, when the geom_sf() is added inside the above it
  ##   throws an error. Temporary solution until fixed:
  p <- p + geom_sf(aes(geometry = geometry, fill = Stat),
                   colour = NA)
  p <- p + geom_sf(aes(geometry = geom_b),
                   fill = NA, colour = signif_col, linewidth = 0.25)

  ## Return plot
  return(p)
}


#'
#' @param clust_col A named character vector of colours to be used for
#' plotting. The names must match the cluster names. Look at details below for
#' more information about what these names should be.
plotSA_localClust <- function(m_sfe,
                              sample_id,
                              feature,
                              statistic = c("moran", "geary", "getis"),
                              test = c("z-score", "permutation"),
                              pVal = 0.05,
                              type = c("spot", "hex"),
                              clust_col = NULL) {
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

  ## Get feature suffixes
  feature <- .int_getSAFeatSuffix(feature = feature, statistic = statistic)

  ## Get the `localResults` name for the specific statistic
  localRes <- .int_getLocalResName(statistic = statistic, test = test)

  ## Check if `localGearyC` is being used. If yes, set `pVal` to NULL
  pVal <- .int_checkSACompatible(localRes = localRes, pVal = pVal)

  ## Prepare the data to plot
  data <- .int_prepareSAPlotData(sfe = sfe,
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

  ## Create the ggplot object
  p <- ggplot2::ggplot(data = data) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(fill = fill)

  ## Add geom_sf
  ## - for an unknown reason, when the geom_sf() is added inside the above it
  ##   throws an error. Temporary solution until fixed:
  p <- p + geom_sf(aes(geometry = geometry, fill = Clust),
                   colour = NA)

  ## Return plot
  return(p)
}


# ---------------------------------------------------------------------------- #
#  ######### INTERNAL FUNCTIONS ASSOCIATED WITH SA PLOTS  (C, G, I) #########
# ---------------------------------------------------------------------------- #
.int_getSAFeatSuffix <- function(feature, statistic) {
  if (statistic == "moran") {
    features <- paste0(feature, c(".Ii", ".IiFDR", ".IiClust"))
  } else if (statistic == "geary") {
    features <- paste0(feature, c(".Ci", ".CiFDR", ".CiClust"))
  } else if (statistic == "getis") {
    features <- paste0(feature, c(".Gi", ".GiFDR", ".GiClust"))
  }

  return(features)
}


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
#' @return
#' This function returns a value for the `pVal` argument.
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


#'
#' @importFrom dplyr case_when
#' @importFrom SpatialFeatureExperiment localResult
#'
.int_prepareSAPlotData <- function(sfe,
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
                                      feature = feature[1]),
                   p.value = localResult(x = sfe,
                                         type = localRes,
                                         feature = feature[2]),
                   geometry = colGeometries(sfe)[[geoms]][["geometry"]]) %>%
          mutate(signif = as.factor(case_when(p.value < pVal ~ "Signif.",
                                              p.value >= pVal ~ "Not Signif.")),
                 geom_b = if_else(signif == "Signif.", .data$geometry, NA))

  if (plot == "cluster") {
    dt <- dt %>%
      mutate(Clust = localResult(x = sfe,
                           type = localRes,
                           feature = feature[3]),
             Clust = as.factor(if_else(signif == "Signif",
                                       .data$Clust,
                                       "Not Signif.")))

  }
  return(dt)
}


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


.int_getSAClustColours <- function(clust_col, clusters) {
  cols_list <- list("Low" = "#0066ff",
                    "High" = "#ff0000",
                    "low-Low" = "#0066ff",
                    "High-Low" = "#ffb3b3",
                    "Low-High" = "#99c2ff",
                    "High-High" = "#ff0000",
                    "Not Signif." = "#E0E0E0")

  if (!is.null(clust_col)) {
    cols_list[names(clust_col)] <- clust_col
  }

  out <- cols_list[clusters]
  out <- unlist(out)

  return(out)

}
