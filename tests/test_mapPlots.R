m_sfe =   msfe
sample_id = sampleNames[1]
type = "hex"
genes = c("ENSG00000254709", "ENSG00000077942")
assay = "logcounts"
minmax = c(0, Inf)
res = "lowres"
fill_args = list(option = "viridis",
                 na.value = "grey")
alpha = 0.3

## Check SFE or MSFE?
sfeT <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

if (type == "hex") {
  stopifnot("spotHex" %in% names(colGeometries(sfeT)))
  type <- "spotHex"
} else if (type == "spot") {
  type <- "spotPoly"
}

## Fetch image if needed
if (!is.null(res)) {
  ## Fetch image data and transform to raster
  image <- .int_getImgDt(sfe = sfeT, sample_id = sample_id, image_id = res)
  ## Get capture area limits
  limits_list <- .int_getImgLims(sfe = sfeT)
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
# genes_df <- assay(sfeT, assay) %>%
#   t() %>%
#   as.data.frame() %>%
#   dplyr::select(tidyr::all_of(genes)) %>%
#   tibble::rownames_to_column(var = "rownames")
genes_df <- assay(sfeT, assay) %>%
  t() %>%
  as.data.frame()

genes_lst <- lapply(genes, FUN = function(x) {
  genes_df %>%
    dplyr::select(tidyr::all_of(x)) %>%
    tibble::rownames_to_column(var = "rownames") %>%
    filter(.data[[x]] > minmax[1] & .data[[x]] < minmax[2])
})

geoms <- data.frame(geometry = colGeometry(sfeT, type)) %>%
  tibble::rownames_to_column(var = "rownames")

geoms_lst <- list()
geoms_lst[[1]] <- geoms

data_lst <- c(geoms_lst, genes_lst)

data <- purrr::reduce(data_lst, dplyr::left_join, by = 'rownames') %>%
  tibble::column_to_rownames(var = "rownames") %>%
  tidyr::pivot_longer(cols = -"geometry",
                      names_to = "gene",
                      values_to = "expression")

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


rm(m_sfe, sfeT, genes_df, fill_args, fill, assay, image, limits_list, type,
   sample_id, genes, assay, minmax, res, alpha, geoms, geoms_lst)
