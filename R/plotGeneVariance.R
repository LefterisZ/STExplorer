plotGeneVariance <- function(dec, hvgs, sample_id = TRUE, ...) {
  ## Check arguments
  # stopifnot(is(msfe, "SpatialFeatureExperiment"))

  ## Select samples
  ids <- .int_getMSFEsmplID(list = dec, sample_id = sample_id)

  ## Get the curve fit and prepare the data frame
  fit_list <- lapply(ids, .int_getFit, dec = dec, hvgs = hvgs)
  fit_df <- rlist::list.rbind(fit_list)

  ## Plot the visualisation
  ggplot(data = fit_df,
         aes(x = mean, y = var, colour = topHVGs)) +
    geom_point() +
    geom_line(aes(y = trend), colour = "dodgerblue", linewidth = 1.5) +
    scale_colour_manual(values = c("black", "red")) +
    facet_wrap(~ sID) +
    labs(x = "Mean of log-expression",
         y = "Variance of log-expression",
         colour = "Top HVGs") +
    theme_classic()
}


.int_getFit <- function(id, dec, hvgs) {
  fit <- metadata(dec[[id]])
  fit_df <- data.frame(mean = fit$mean,
                       var = fit$var,
                       trend = fit$trend(fit$mean),
                       sID = id)

  fit_df <- fit_df %>%
    tibble::rownames_to_column(var = "row.names") %>%
    dplyr::mutate(topHVGs = ifelse(row.names %in% hvgs[[id]], TRUE, FALSE)) # %>%
    # tibble::column_to_rownames("row.names")

  return(fit_df)
}
