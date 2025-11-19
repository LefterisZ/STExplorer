#' Get MSigDB Data
#'
#' This function retrieves data from MSigDB for a specified species.
#'
#' @param species A character string specifying the target species
#' (e.g., "Homo sapiens").
#' @param ... Other arguments passed onto the \code{msigdbr} function.
#'
#' @return A data table with the MSigDB data.
#'
#' @importFrom msigdbr msigdbr
#' @importFrom data.table data.table
#'
#' @examples
#' # Download data from MSigDB
#' msigdb <- getMSigDBData("Homo sapiens")
#'
#' @export
getMSigDBData <- function(species, ...) {
  ## Download data from MSigDB
  msig_data <- data.table::data.table(msigdbr::msigdbr(species = species, ...))

  return(msig_data)
}

#' Look at Gene Set categories and their subcategories
#'
#' This function is a wrapper around the "msigdbr_collections" function that
#' prints the possible combinations of categories and subcategories. It is
#' useful because if a category is matched with a subcategory that it shouldn't
#' then GSEA will return no results.
#'
#' @importFrom msigdbr msigdbr_collections
#'
#' @examples
#' # Check available collections
#' viewCollections()
#'
#' @export
viewCollections <- function(){
  ## Check available collections
  print(n = 30, msigdbr::msigdbr_collections())
}


#' Get Term-to-Gene Data
#'
#' This function retrieves and prepares a term-to-gene data frame based on user
#' choices.
#'
#' @param user_data A user-supplied term-to-gene data frame (optional).
#' @param msig_data A term-to-gene data frame as generated using the
#' `getMSigDBData` function.
#' @param cat Filter for "gs_cat" (e.g., "C2", "C3").
#' @param subcat Filter for "gs_subcat" (e.g., "GO:BP", "CGP").
#'
#' @return A term-to-gene data frame with columns "term" and "gene."
#'
#' @importFrom data.table data.table
#' @importFrom data.table .SD
#'
#' @examples
#' # Get MSigDB data
#' data <- getMSigDBData("Homo sapiens")
#'
#' # View categories and subcategories
#' viewCollections()
#'
#' # Get term-to-gene data.table
#' t2g <- getTerm2Gene(msig_data = data, cat = C3, subcat = "")
#'
#' # Note in the above that if a category has no subcategory then we put "" in
#' # the 'subcat' argument.
#'
#' #' @seealso \code{\link[clusterProfiler]{GSEA}}
#'
#' @export
getTerm2Gene <- function(user_data = NULL,
                         msig_data = NULL,
                         cat = NULL,
                         subcat = NULL) {
  if (!is.null(user_data)) {
    ## User provided their own data frame
    ## Perform checks to ensure it has "term" and "gene" columns
    if (!all(c("term", "gene") %in% colnames(user_data))) {
      stop("User-supplied data frame must have 'term' and 'gene' columns.")
    }
    return(data.table(user_data))

  } else if (!is.null(msig_data)) {
    ## User chooses to use MSigDB data
    msig_data <- data.table::data.table(msig_data)
    ## Perform subsetting based on "gs_cat" and "gs_subcat"
    if (!is.null(cat) && !is.null(subcat)) {
      msig_data <- msig_data[msig_data$gs_collection %in% cat &
                               msig_data$gs_subcollection %in% subcat]
    }

    ## Perform subsetting to select "gs_name" and the appropriate
    ## "ensembl_gene" column
    gene_column <- grep("_ensembl_gene", colnames(msig_data), value = TRUE)
    result <- msig_data[, .(term = .SD[["gs_name"]],
                            gene = .SD[[gene_column]]),
                        .SDcols = c("gs_name", gene_column)]
    return(result)

  } else {
    stop("You must provide user data or specify an MSigDB data frame.",
         "\nFor the later, we recommend using the 'getMSigDBData' function")
  }
}



#' Run Gene Set Enrichment Analysis (GSEA) per location using GWPCA
#'
#' This function performs GSEA on multiple cores to analyze the enrichment of
#' gene sets in different spatial locations based on Gene Set Enrichment
#' Analysis. It takes the output of a GWPCA analysis, extracts the necessary
#' data for each spatial location, and applies GSEA. The results are combined
#' into a final dataframe, and additional filtering and processing are applied.
#'
#' @param gwpca A GWPCA object containing the results of the spatial
#' transcriptomics analysis.
#' @param pc The principal component (PC) index to be used for GSEA.
#' @param genes_no The minimum number of genes from the gwpca results that has
#' to be present in a gene set to be considered for enrichment.
#' @param NES The minimum Normalized Enrichment Score (NES) for considering a
#' gene set as enriched.
#' @param minGSSize The minimum gene set size to be considered in GSEA.
#' @param pvalueCutoff The p-value cutoff for identifying significant gene sets
#' in GSEA.
#' @param TERM2GENE The term-to-gene mapping for gene sets.
#' @param pAdjustMethod The method for multiple testing correction in GSEA
#' (default is "fdr").
#' @param scoreType The GSEA scoring type (default is "pos"). Possible options
#'                  are ("std", "pos", "neg"). If the "std" option is selectedthe
#'                  enrichment score is computed as in the original GSEA. The
#'                  "pos" and "neg" score types are intended to be used for
#'                  one-tailed tests (i.e. when one is interested only in
#'                  positive ("pos") or negateive ("neg") enrichment). Look at
#'                  details for more info.
#' @param nPermSimple The number of permutations for simple GSEA.
#' @param mc.cores The number of cores to use for parallel processing. If not
#' provided and type is "FORK", it defaults to availableCores() - 1.
#' @param type The type of cluster to be created. Default is "FORK". Other
#' options may include "PSOCK" or others supported by parallel::makeCluster.
#' @param regex A regular expression pattern for cluster identification.
#' @param ... Additional arguments to be passed to internal functions.
#'
#' @details
#' Detailed explanation of the function, including any relevant details about the
#' algorithm used, data processing steps, and result interpretation.
#'
#' For the score type the default is "pos" because here we use PCA loadings
#' instead of differential gene expression. The ± signs in a PCA loading are
#' somewhat arbitrary and depend on the initial rotation used in PCA
#' calculations. For this reason, we suggest the "pos" approach where we
#' provide to GSEA the absolute loadings and then essentially we look only for
#' positive enrichment. You are free to select the "std" method. In this case,
#' the original loadings are passed into GSEA.
#'
#' @importFrom stringr str_count
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @importFrom dplyr bind_rows group_by filter arrange slice_max mutate
#' @importFrom tibble column_to_rownames
#'
#' @return A tibble containing the GSEA results for each spatial location.
#'
#' @keywords spatial transcriptomics gene set enrichment analysis GSEA GWPCA
#'
#' @family Spatial Transcriptomics Analysis
#'
#' @aliases gwpca_FunctionalClustering
#'
#' @rdname gwpca_FunctionalClustering
#'
#' @examples
#' ## Example usage:
#' # gwpca_result <- gwpcaSTE()  # run the gwpcaSTE function
#'
#' ## Functional clustering on Principal Component 1 (PC1)
#' # gsea_result <- gwpca_FunctionalClustering(gwpca_result, pc = 1)
#'
#' @seealso \code{\link{gwpcaSTE}}, \code{\link[clusterProfiler]{GSEA}}
#'
#' @export
gwpca_FunctionalClustering <- function(gwpca,
                                       pc,
                                       genes_no = 1,
                                       NES = 1.5,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 4,
                                       type = "FORK",
                                       regex = "",
                                       ...) {
  ## Get input dataset from gwpca
  inputGSEA <- gwpca$loadings[,,pc]

  ## Run GSEA for all locations
  # list <- mclapply(X = 1:nrow(inputGSEA),
  #                  FUN = .int_runGSEA,
  #                  inputGSEA = inputGSEA,
  #                  minGSSize = minGSSize,
  #                  pvalueCutoff = pvalueCutoff,
  #                  TERM2GENE = TERM2GENE,
  #                  pAdjustMethod = pAdjustMethod,
  #                  scoreType = scoreType,
  #                  nPermSimple = nPermSimple,
  #                  mc.cores = mc.cores,
  #                  ...)

  my.cl <- makeClusterGWPCA(spec = mc.cores, type = type)
  doParallel::registerDoParallel(my.cl)

  list <- foreach::foreach(X = 1:nrow(inputGSEA)) %dopar% {
    STExplorer:::.int_runGSEA(X = X,
                              inputGSEA = inputGSEA,
                              minGSSize = minGSSize,
                              pvalueCutoff = pvalueCutoff,
                              TERM2GENE = TERM2GENE,
                              pAdjustMethod = pAdjustMethod,
                              scoreType = scoreType,
                              nPermSimple = nPermSimple,
                              ...)
  }

  parallel::stopCluster(my.cl)

  ## Make the dataframe
  gsea <- bind_rows(list)

  genesNo <- genes_no
  NES_thresh <- NES

  ## Check if the gsea data frame is populated. If it is not, the next code
  ## chunk will error with:
  # `dplyr::mutate()`:
  # ℹ In argument: `genes_no = str_count(.data$core_enrichment, "ENSG")`.
  # Caused by error in `.data$core_enrichment`:
  #   ! Column `core_enrichment` not found in `.data`.
  if (sum(is.na(gsea$ID)) == nrow(gsea)) {
    stop("The GSEA returned no results in any location of the sample.\n",
         "Please check that gene ids from GWPCA are contained in the \n",
         "TERM2GENE list you provide. If this is not the source, consider ",
         "further lowering the pValue threshold.")
  }

  gsea <-  gsea %>%
    dplyr::mutate(genes_no = str_count(.data$core_enrichment,"ENS")) %>%
    dplyr::group_by(.data$Location) %>%
    dplyr::filter(genes_no > genesNo) %>%
    dplyr::filter(abs(NES) > NES_thresh, .by_group = TRUE) %>%
    # dplyr::filter(qvalue < 0.3, .by_group = TRUE) %>%
    # dplyr::arrange(rank, .by_group = TRUE) %>%
    # dplyr::slice_min(rank, n = 1, with_ties = FALSE) %>%
    dplyr::arrange(NES, .by_group = TRUE) %>%
    dplyr::slice_max(NES, n = 1, with_ties = FALSE) %>%
    tibble::column_to_rownames(var = "Location")

  regex <- regex
  gsea_map <- merge(gsea,
                    gwpca$geometry,
                    by = "row.names",
                    all.y = TRUE) %>%
    dplyr::mutate(cluster = gsub(regex, "", .data$ID)) %>%
    tibble::column_to_rownames("Row.names")

  return(gsea_map)
}


# ---------------------------------------------------------------------------- #
#  ################# INTERNAL FUNCTIONS ASSOCIATED WITH GWR #################
# ---------------------------------------------------------------------------- #
#' Internal Function: Run Gene Set Enrichment Analysis (GSEA) for a Single
#' Location
#'
#' This function performs GSEA for a single spatial location using a given set
#' of parameters. It is an internal function used within the main analysis
#' function and is not intended for direct use.
#'
#' @param X The index of the spatial location for which GSEA is performed.
#' @param inputGSEA The input gene expression data for the specified spatial
#' location.
#' @param minGSSize The minimum gene set size to be considered in GSEA.
#' @param pvalueCutoff The p-value cutoff for identifying significant gene sets
#' in GSEA.
#' @param TERM2GENE The term-to-gene mapping for gene sets.
#' @param pAdjustMethod The method for multiple testing correction in GSEA.
#' @param verbose A logical value indicating whether to print verbose output
#' during GSEA.
#' @param scoreType The GSEA scoring type.
#' @param nPermSimple The number of permutations for simple GSEA.
#' @param ... Additional arguments to be passed to the GSEA function.
#'
#' @importFrom clusterProfiler GSEA
#'
#' @return A dataframe containing the results of GSEA for the specified spatial
#' location.
#'
#' @examples
#' # This function is an internal function used within the main function.
#' # Refer to the main analysis function for example usage.
#'
#' @keywords internal gene set enrichment analysis GSEA spatial transcriptomics
#' @family Spatial Transcriptomics Analysis
#'
#' @aliases .int_runGSEA
#'
#' @rdname dot-int_runGSEA
#'
#' @details
#' This function is used internally to perform Gene Set Enrichment Analysis
#' (GSEA) for a single spatial location based on the provided gene expression
#' data and parameters. The results are returned as a dataframe.
#'
#' @seealso
#' \code{\link{gwpca_FunctionalClustering}}, \code{\link[clusterProfiler]{GSEA}}
#'
.int_runGSEA <- function(X,
                         inputGSEA,
                         minGSSize,
                         pvalueCutoff,
                         TERM2GENE,
                         pAdjustMethod,
                         verbose,
                         scoreType,
                         nPermSimple,
                         ...) {
  # Generate geneList
  geneList <- .int_generateGeneList(inputGSEA = inputGSEA,
                                    X = X,
                                    scoreType = scoreType)

  # Run GSEA
  gsea <- clusterProfiler::GSEA(geneList,
                                minGSSize = minGSSize,
                                pvalueCutoff = pvalueCutoff,
                                TERM2GENE = TERM2GENE,
                                pAdjustMethod = pAdjustMethod,
                                verbose = FALSE,
                                scoreType = scoreType,
                                nPermSimple = nPermSimple,
                                ...)

  if (dim(gsea@result)[1] == 0) {
    na <- c(rep(NA, dim(gsea@result)[2])) %>%
      t() %>%
      as.data.frame()
    colnames(na) <- colnames(gsea@result)
    gsea@result <- rbind(gsea@result, na)
  }
  # Add location barcode
  gsea@result$Location <- rownames(inputGSEA)[X]

  return(gsea@result)
}


#' Internal Function: Generate Gene List Based on Score Type
#'
#' This function generates a gene list based on the specified score type from
#' the input gene expression data for a single spatial location. It is an
#' internal function used within the main analysis function and is not intended
#' for direct use.
#'
#' @param inputGSEA The input gene expression data for the specified spatial
#' location.
#' @param X The index of the spatial location for which the gene list is
#' generated.
#' @param scoreType The type of score used for sorting genes. Valid options
#' are "std" for standard GSEA and "pos" for positive enrichment.
#'
#' @return A vector containing the sorted gene list.
#'
#' @aliases .int_generateGeneList
#'
#' @rdname dot-int_generateGeneList
#'
#' @details
#' This function internally sorts genes in descending order based on the
#' specified score type ("std" or "pos"). For "std", genes are sorted based on
#' their original loading scores, while for "pos", genes are sorted based on
#' the absolute value of their expression values.
#'
#' @seealso
#' \code{\link{.int_sortGenes}}
#'
#' @keywords internal gene expression data gene list sorting score type
#' @family Spatial Transcriptomics Analysis
#'
#' @examples
#' # This function is an internal function used within the main function.
#' # Refer to the main analysis function for example usage.
#'
.int_generateGeneList <- function(inputGSEA, X, scoreType) {
  if (scoreType == "std") {
    geneList <- .int_sortGenes(inputGSEA[X,])
  } else if (scoreType == "pos") {
    geneList <- .int_sortGenes(abs(inputGSEA[X,]))
  }
  return(geneList)
}


#' Internal Function: Sort Genes in Descending Order
#'
#' This function sorts genes in descending order based on their expression
#' values. It is an internal function used within other functions and is not
#' intended for direct use.
#'
#' @param data The gene expression data to be sorted.
#'
#' @return A vector containing the sorted gene expression values.
#'
#' @aliases .int_sortGenes
#'
#' @rdname dot-int_sortGenes
#'
#' @details
#' This function sorts genes in descending order based on their expression
#' values. It is used internally within other functions to sort genes before
#' further processing.
#'
#' @seealso
#' \code{\link{.int_generateGeneList}}
#'
#' @keywords internal gene expression data gene list sorting
#'
#' @examples
#' # This function is an internal function used within other functions.
#' # Refer to the documentation of the calling function for example usage.
#'
.int_sortGenes <- function(data) {
  return(data[order(data, decreasing = TRUE)])
}



# additional helpers - not exported
keepSize_t2g <- function(t2g_df, min_size = 15, max_size = 500) {
  stopifnot(all(c("term", "gene") %in% names(t2g_df)))

  t2g_df %>%
    distinct(term, gene) %>%
    add_count(term, name = "term_size") %>%
    filter(term_size >= min_size, term_size <= max_size) %>%
    select(-term_size)
}


# t2g: data.frame with columns term, gene (case-insensitive accepted)
# keep_patterns / drop_patterns: character vectors of regexes
# allowed_terms: optional exact term names to keep (vector of strings)
# size filter applied AFTER pattern filtering
filter_t2g <- function(t2g,
                       keep_patterns = character(),
                       drop_patterns = character(),
                       allowed_terms = NULL,
                       min_size = 10,
                       max_size = Inf) {
  # normalise columns
  nm <- tolower(names(t2g))
  names(t2g) <- nm
  if (!all(c("term","gene") %in% names(t2g))) {
    # common fallbacks
    col_map <- c(gs_name = "term", term_name = "term", pathway = "term", set = "term", name = "term",
                 gene_symbol = "gene", genesymbol = "gene", symbol = "gene", geneid = "gene",
                 entrez = "gene", gene_id = "gene")
    for (k in names(col_map)) if (k %in% names(t2g)) names(t2g)[names(t2g) == k] <- col_map[[k]]
    stopifnot(all(c("term","gene") %in% names(t2g)))
  }

  t2g <- t2g %>% distinct(term, gene)

  # If an explicit allow-list is provided, use it (exact match)
  if (!is.null(allowed_terms)) {
    t2g <- t2g %>% semi_join(tibble(term = allowed_terms), by = "term")
  } else {
    # otherwise apply regex keep/drop
    any_match <- function(x, patterns, default = TRUE) {
      if (length(patterns) == 0) return(rep(default, length(x)))
      str_detect(x, regex(paste(patterns, collapse = "|"), ignore_case = TRUE))
    }
    tl <- tolower(t2g$term)
    keep_flag <- any_match(tl, keep_patterns, default = TRUE)  # TRUE if no keep patterns
    drop_flag <- any_match(tl, drop_patterns, default = FALSE)
    t2g <- t2g[keep_flag & !drop_flag, , drop = FALSE]
  }

  # size limits
  sizes <- t2g %>% count(term, name = "set_size")
  keep_terms <- sizes %>% filter(set_size >= min_size, set_size <= max_size) %>% pull(term)
  t2g %>% semi_join(tibble(term = keep_terms), by = "term")
}


pretty_term_labels <- function(terms,
                               remove_prefix = TRUE,
                               title_case = TRUE,
                               keep_go_id = FALSE) {
  out <- terms

  if (remove_prefix) {
    # strip common MSigDB namespaces / collection tags at the start
    out <- out %>%
      str_replace("^C[0-9]+[:_]", "") %>%                 # C2_, C5_, etc (rare in names)
      str_replace("^(GO|GOBP|GOMF|GOCC)[:_]", "") %>%
      str_replace("^REACTOME[:_]", "") %>%
      str_replace("^HALLMARK[:_]", "") %>%
      str_replace("^CP[:_]", "") %>%
      str_replace("^WP[:_]", "") %>%
      str_replace("^KEGG[:_]", "") %>%
      str_replace("^M[0-9]+[:_]", "") %>%                 # M8_, etc
      str_replace("^MOUSE[:_]", "")
  }

  # Optionally keep GO ID at the end as "(GO:0001234)" if present in the original
  if (keep_go_id) {
    go_id <- str_extract(terms, "GO:\\d{7}")
    out <- ifelse(!is.na(go_id), paste0(out, " (", go_id, ")"), out)
  }

  # underscores -> spaces; squeeze whitespace
  out <- out %>%
    str_replace_all("_", " ") %>%
    str_squish()

  if (title_case) {
    out <- str_to_title(out)
    # preserve common acronyms
    fix <- function(x, patt) str_replace_all(x, patt)
    out <- fix(out, c(
      "\\bNmda\\b" = "NMDA",
      "\\bAmpa\\b" = "AMPA",
      "\\bGaba\\b" = "GABA",
      "\\bMtorc1\\b" = "mTORC1",
      "\\bBbb\\b" = "BBB",
      "\\bDna\\b" = "DNA",
      "\\bRna\\b" = "RNA"
    ))
  }

  out
}
