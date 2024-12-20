#' Find Leading genes
#'
#' @name gwpca_LeadingGene
#'
#' @description
#' This function identifies the leading genes in a location for one or multiple
#' principal components.
#'
#' @param gwpca A list of class \code{gwpca} containing the geographically
#' weighted principal component analysis results.
#'
#' @param m_sfe A \code{SpatialFeatureExperiment} or a
#' \code{MetaSpatialFeatureExperiment} object with the spatial
#' feature experiment data.
#'
#' @param sample_id The sample ID (or the SFE object name inside a Meta-SFE
#' object) for which the GWPCA was run.
#'
#' @param pc_nos A vector containing the principal component numbers for which
#' you want to find the leading genes.
#'
#' @param genes_n An integer indicating how many genes you want to include in
#' the top leading genes. This argument is used only when the type parameter is
#' set to 'multi'. The default value is 5.
#'
#' @param type A character indicating the type of leading genes to be
#' identified. Possible values are 'single' for a single leading gene per
#' location or 'multi' for multiple leading genes per location.
#'
#' @param method A character indicating the method to be used for grouping the
#' spots together. Possible values are "membership" or "order".
#'
#' @param sort A character indicating whether the absolute or the original
#' leading score will be used for sorting. Possible values are "abs" (for
#' absolute) or "original".
#'
#' @param names A character indicating the output format for gene identifiers.
#' Possible values are "id" for ENSG IDs or "gene_names" for gene names.
#'
#' @details
#' The function calculates the leading genes in a location based on the
#' specified principal component numbers. It returns the results as a list of
#' leading genes for each principal component. If the type parameter is set to
#' 'single', it returns a single leading gene for each location and principal
#' component. If the type parameter is set to 'multi', it returns multiple
#' leading genes (specified by genes_n) for each location and principal
#' component.
#'
#' @returns A modified \code{gwpca} object with additional information on the
#' leading genes. If the type parameter is set to 'single', the modified object
#' contains a new slot leadingGeneSingle with the leading genes for each
#' principal component. If the type parameter is set to 'multi', the modified
#' object contains a new slot leadingGeneMulti with the top leading genes
#' (specified by genes_n) for each principal component. The leading genes are
#' represented as a data frame with columns corresponding to principal component
#' numbers and rows corresponding to locations. The data frame also includes the
#' geometry information of the locations.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @importFrom SpatialFeatureExperiment colGeometry
#' @importFrom dplyr select starts_with matches
#'
#' @examples
#' # Load Geographically Weighted principal component analysis data
#' data(gwpca)
#' data(sfe)
#'
#' #Identify single leading gene per location for PC1 and PC2
#' gwpca <- gwpca_LeadingGene(gwpca = gwpca, sfe = sfe, pc_nos = c(1, 2),
#' type = "single", names = "id")
#'
#' # Identify top 3 leading genes per location for PC1 and PC2
#' gwpca <- gwpca_LeadingGene(gwpca = gwpca, sfe = sfe, pc_nos = c(1, 2),
#' type = "multi", genes_n = 3, method = "membership", names = "gene_names")
#'
#' @export
gwpca_LeadingGene <- function(gwpca,
                              m_sfe,
                              sample_id = NULL,
                              pc_nos,
                              genes_n = 5,
                              type = c("single", "multi"),
                              method = c("membership", "order"),
                              sort = c("abs", "original"),
                              names = c("id", "gene_names")) {
  ## Check valid argument inputs
  type <- match.arg(type)
  method <- match.arg(method)
  sort <- match.arg(sort)
  names <- match.arg(names)

  ## SFE or metaSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Find the leading genes for multiple principal components
  if (type == "single") {
    leading_genes <- lapply(pc_nos, function(pc_no) {
      int_LeadingGene_single(gwpca, pc_no, sfe, sort)
    })
  } else if (type == "multi") {
    leading_genes <- lapply(pc_nos, function(pc_no) {
      int_LeadingGene_multi(gwpca, pc_no, genes_n, sfe, method, sort)
    })
  }

  ## Concatenate to a data frame
  result <- do.call(cbind.data.frame, leading_genes)
  if (names == "id") {
    colnames(result)[colnames(result) == "id"] <- paste0("id", pc_nos)
  } else {
    colnames(result)[colnames(result) == "gene_names"] <- paste0("gene_names",
                                                                 pc_nos)
  }

  result$geometry <- colGeometry(sfe, "spotHex")

  ## Select ENSG IDs or gene names to be in the final output
  result <- dplyr::select(result, starts_with(names) | matches("geometry"))
  colnames(result)[-ncol(result)] <- paste0("PC", pc_nos)

  ## Add colnames
  if (type == "single") {
    gwpca$leadingGeneSingle <- result
  } else if (type == "multi") {
    gwpca$leadingGeneMulti <- result
  }

  return(gwpca)
}

# ----------------------------NON-EXPORTED FUNCTIONS-------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Find the single leading gene
#'
#' @name int_LeadingGene_single
#'
#' @description
#' A function to find the SINGLE leading gene per principal component in a
#' location.
#'
#' @keywords internal
#'
#' @param gwpca a list of class \code{gwpca}.
#' @param pc_no a numeric value of the principal component number for which
#' you want to find the leading genes.
#' @param sfe a \code{SpatialFeatureExperiment} object.
#' @param sort A character indicating whether the absolute or the original
#' leading score will be used for sorting. Possible values are "abs" (for
#' absolute) or "original"
#'
#' @importFrom SpatialFeatureExperiment rowData
#'
int_LeadingGene_single <- function(gwpca, pc_no, sfe, sort) {
    pc_name <- paste0("PC", pc_no)

    local_loadings <- round(gwpca$loadings[, , pc_no], 4)

    if (sort == 'abs') {
      local_loadings <- abs(local_loadings)
    }

    leading_gene_indices <- max.col(local_loadings)
    leading_gene_ids <- colnames(local_loadings)[leading_gene_indices]

    if (is.null(rowData(sfe)$id)) {
        gene_name_indices <- match(leading_gene_ids, rownames(sfe))
    } else {
        gene_name_indices <- match(leading_gene_ids, rowData(sfe)$id)
    }
    leading_gene_names <- rowData(sfe)$gene_name[gene_name_indices]

    result <- data.frame(id = leading_gene_ids,
                         gene_names = leading_gene_names,
                         stringsAsFactors = FALSE)

    leading_gene_counts <- table(result$gene_names)
    cat(length(unique(leading_gene_ids)),
        " leading genes found for ", pc_name)
    cat("\nThe leading genes in ", pc_name, " are:")
    print(leading_gene_counts)

    return(result)
}
# ---------------------------------------------------------------------------- #
#' Find the top-k leading genes
#'
#' @name int_LeadingGene_multi
#'
#' @description
#' A function to find the top k leading genes per principal components in a
#' location.
#'
#' @keywords internal
#'
#' @param gwpca a list of class \code{gwpca}.
#' @param pc_no a numeric value for the principal component numbers for which
#' you want to find the leading genes.
#' @param genes_n an integer indicating how many genes you want to be included.
#' @param sfe a \code{SpatialFeatureExperiment} object.
#' @param method takes values either "membership" or "order". Is the method to
#' be used for grouping the spots together.
#' @param sort A character indicating whether the absolute or the original
#' leading score will be used for sorting. Possible values are "abs" (for
#' absolute) or "original"
#'
#' @importFrom SpatialFeatureExperiment rowData
#'
int_LeadingGene_multi <- function(gwpca, pc_no, genes_n, sfe, method, sort) {
    pc_names <- paste0("PC", pc_no)

    local_loadings <- round(gwpca$loadings[, , pc_no], 4)

    if (sort == 'abs') {
      local_loadings <- abs(local_loadings)
    }

    ## Order loading scores in decreasing order
    ordered_indices <- apply(local_loadings, 1, order, decreasing = TRUE)

    ## Get the top leading genes by membership or order
    top_genes <- apply(ordered_indices[1:genes_n, ], 2, function(indices) {
        gene_ids <- colnames(local_loadings)[indices]
        if (is.null(rowData(sfe)$id)) {
            gene_name_indices <- match(gene_ids, rownames(sfe))
        } else {
            gene_name_indices <- match(gene_ids, rowData(sfe)$id)
        }
        gene_names <- rowData(sfe)$gene_name[gene_name_indices]
        if (method == "membership") {
            gene_ids <- paste(sort(gene_ids), collapse = ";")
            gene_names <- paste(sort(gene_names), collapse = ";")
        } else {
            gene_ids <- paste(gene_ids, collapse = ";")
            gene_names <- paste(gene_names, collapse = ";")
        }
        data.frame(id = gene_ids,
                   gene_names = gene_names,
                   stringsAsFactors = FALSE)
    })

    leading_genes <- do.call(rbind, top_genes)


    unique_groups <- unique(leading_genes$gene_names)
    cat("The number of individual leading genes groups found for", pc_names,
        "is:", length(unique_groups), "\nThese groups are:")
    if (length(unique_groups) < 15) {
        print(table(leading_genes$gene_names))
    } else {
        cat(" Too many to print them!\n")
    }

    ## Return
    return(leading_genes)
}
