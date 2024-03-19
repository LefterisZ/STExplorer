# pc = 1,
genes_no = 2
NES = 1.5
# minGSSize = 5,
# pvalueCutoff = 0.25,
# TERM2GENE,
# pAdjustMethod = "fdr",
# scoreType = "std",
# nPermSimple = 10000,
# mc.cores = 4,
# regex = ""

## Get input dataset from gwpca
inputGSEA <- gwpca$loadings[,,1]

## Run GSEA for all locations
list <- parallel::mclapply(X = 1:nrow(inputGSEA),
                           FUN = STExplorer:::.int_runGSEA,
                           inputGSEA = inputGSEA,
                           minGSSize = 5,
                           pvalueCutoff = 0.25,
                           TERM2GENE = t2g,
                           pAdjustMethod = "fdr",
                           scoreType = "std",
                           nPermSimple = 10000,
                           mc.cores = 8)

## Make the dataframe
gsea <- dplyr::bind_rows(list)

genesNo <- genes_no
NES_thresh <- NES

gsea <-  gsea %>%
  dplyr::mutate(genes_no = stringr::str_count(.data$core_enrichment,"ENSG")) %>%
  dplyr::group_by(.data$Location) %>%
  dplyr::filter(genes_no > genesNo) %>%
  dplyr::filter(abs(NES) > NES_thresh, .by_group = TRUE) %>%
  # dplyr::filter(qvalue < 0.3, .by_group = TRUE) %>%
  # dplyr::arrange(rank, .by_group = TRUE) %>%
  # dplyr::slice_min(rank, n = 1, with_ties = FALSE) %>%
  dplyr::arrange(NES, .by_group = TRUE) %>%
  dplyr::slice_max(NES, n = 1, with_ties = FALSE) %>%
  tibble::column_to_rownames(var = "Location")

regex <- ""
gsea_map <- merge(gsea,
                  gwpca$geometry,
                  by = "row.names",
                  all.y = TRUE) %>%
  dplyr::mutate(cluster = gsub(regex, "", .data$ID)) %>%
  tibble::column_to_rownames("Row.names")


# rm(inputGSEA, list, gsea, regex, genes_no, NES)
