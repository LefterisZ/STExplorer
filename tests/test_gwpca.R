gwpca = gwpca
comps = c(1,3,5)
genes = c("ENSG00000254709", "ENSG00000077942", "ENSG00000107317")
type = "hex"
colours = "viridis"
col_args = list()
gene_names = c("A", "B", "C")

#--------------------------------------------------#
#--------------------------------------------------#
genes_ar <- gwpca$loadings[,genes,comps]

if (!is.null(gene_names)) {
  if (is.character(gene_names)) {
    colnames(genes_ar) <- gene_names
  }
}

genes_lst <- lapply(seq_along(comps), FUN = function(x) {
  genes_ar[,,x] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "rownames") %>%
    dplyr::mutate(pcs = dimnames(genes_ar)[[3]][x])
})

genes_dt <- rlist::list.rbind(genes_lst)

geoms <- gwpca$geometry %>%
  tibble::rownames_to_column(var = "rownames")

out <- dplyr::right_join(geoms, genes_dt, by = 'rownames') %>%
  tidyr::pivot_longer(cols = -c("rownames", "geometry", "pcs"),
                      names_to = "gene",
                      values_to = "leadScore") %>%
  dplyr::mutate(leadScore = abs(leadScore))

data <- out

#--------------------------------------------------#
#--------------------------------------------------#
## Fetch the data
data <- .int_prepDtMapGWPCA(gwpca = gwpca,
                            locs = NULL,
                            comps = comps,
                            genes = genes,
                            type = type,
                            gene_names = gene_names)
if (type == "hex") {
  stopifnot("spotHex" %in% names(colGeometries(sfe)))
  type <- "spotHex"
} else if (type == "spot") {
  type <- "spotPoly"
}

## Set some defaults if not provided

if (DelayedArray::isEmpty(col_args) & colours == "viridis") {
  col_args <- list(option = "viridis")
} else if (DelayedArray::isEmpty(col_args) & colours == "custom") {
  stop("When selecting custom colours, `col_args` must not be empty.")
}

## Generate plot
p <- ggplot(data) +
  ggplot2::geom_sf(aes(geometry = geometry, fill = leadScore)) +
  ggplot2::facet_grid(rows = ggplot2::vars(data$gene),
                      cols = ggplot2::vars(data$pcs))

if (colours == "viridis") {
  p <- p + do.call(scale_fill_viridis_c, c(list(), col_args))
} else if (colours == "custom") {
  p <- p + do.call(scale_fill_gradientn, c(list(), col_args))
}

p +
  ggplot2::labs(fill = "Absolute\nLeading\nscore") +
  ggplot2::coord_sf() +
  ggplot2::theme_void()

rm(comps, genes, type, p, col_args, genes_ar, genes_lst, geoms, geoms_lst, colours, data_lst, out, genes_dt, data, organism, gene_name,  gene_names, mart)
