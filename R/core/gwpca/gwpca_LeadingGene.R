#' @name gwpca_LeadingGene
#' @description
#' A function that identifies the SINGLE leading gene in a location or the k 
#' leading genes in a location.
#' 

# ---------------------------------------------------------------------------- #
#' @name int_LeadingGene_single
#' @description
#' A function to find the SINGLE leading gene per principal component in a 
#' location.
#' 
#' @param gwpca a list of class \code{gwpca}.
#' @param pc.no a numeric value of the principal component number for which 
#' you want to find the leading genes.
#' @param sfe a \code(SpatialFeatureExperiment) object.
#' 
int_LeadingGene_single <- function(gwpca, pc_no, sfe) {
    pc_name <- paste0("PC", pc_no)
    
    local_loadings <- round(gwpca$loadings[, , pc_no], 4)
    
    leading_gene_indices <- max.col(abs(local_loadings))
    leading_gene_ids <- colnames(local_loadings)[leading_gene_indices]
    
    gene_name_indices <- match(leading_gene_ids, rowData(sfe)$id)
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
#' @name int_LeadingGene_multi
#' @description
#' A function to find the top k leading genes per principal components in a 
#' location.
#' 
#' @param gwpca a list of class \code{gwpca}.
#' @param pc.no a numeric value for the principal component numbers for which 
#' you want to find the leading genes.
#' @param genes.n an integer indicating how many genes you want to be included.
#' @param sfe a \code(SpatialFeatureExperiment) object.
#' @param method takes values either "membership" or "order". Is the method to 
#' be used for grouping the spots together.
int_LeadingGene_multi <- function(gwpca, pc_no, genes_n, sfe, method) {
    pc_names <- paste0("PC", pc_no)
    
    local_loadings <- round(gwpca$loadings[, , pc_no], 4)
    
    ## Order loading scores in decreasing order
    ordered_indices <- apply(abs(local_loadings), 1, order, decreasing = TRUE)
    
    ## Get the top leading genes by membership or order
    top_genes <- apply(ordered_indices[1:genes_n, ], 2, function(indices) {
        gene_ids <- colnames(local_loadings)[indices]
        gene_names <- rowData(sfe)$gene_name[match(gene_ids, rowData(sfe)$id)]
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


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' @name gwpca_LeadingGene
#' @description
#' Identify the leading genes in a location for one or multiple principal components.
#' 
#' @param gwpca a list of class \code{gwpca}.
#' @param sfe a \code(SpatialFeatureExperiment) object.
#' @param pc.nos a vector containing the principal component numbers for which 
#' you want to find the leading genes.
#' @param genes.n an integer indicating how many genes you want to be included
#' in the k leading genes. Used only when `type` argument is 'multi'. Defaults 
#' to 5 genes.
#' @param type either 'single' for single leading gene per location or 'multi'
#' for k leading genes per location.
#' @param method takes values either "membership" or "order". Is the method to 
#' be used for grouping the spots together.
#' @param names do you like the output be in ENSG IDs or gene names?
#' 
#' @export

gwpca_LeadingGene <- function(gwpca, 
                              sfe, 
                              pc_nos, 
                              genes_n = 5, 
                              type = c("single", "multi"), 
                              method = c("membership", "order"), 
                              names = c("id", "gene_names")) {
    ## Check valid argument inputs
    type <- match.arg(type)
    method <- match.arg(method)
    names <- match.arg(names)
    
    ## Find the leading genes for multiple principal components
    if (type == "single") {
        leading_genes <- lapply(pc_nos, function(pc_no) {
            int_LeadingGene_single(gwpca, pc_no, sfe)
        })
    } else if (type == "multi") {
        leading_genes <- lapply(pc_nos, function(pc_no) {
            int_LeadingGene_multi(gwpca, pc_no, genes_n, sfe, method)
        })
    }
    
    ## Concatenate to a data frame
    result <- do.call(cbind.data.frame, leading_genes)
    colnames(result)[colnames(result) == "id"] <- paste0("id", pc_nos)
    colnames(result)[colnames(result) == "gene_names"] <- paste0("gene_names", 
                                                                 pc_nos)
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
