#' @name gwpca.topLeading.Gs 
#' 
#' @param gwpca a gwpca output.
#' @param pc.no the Principal Component (PC) to be examined.
#' @param genes.n an integer indicating how many genes you want to be included.
#' @sf.geom the geometry column from an sf object containing a polygons 
#'          geometry column to be used for mapping. It can be a vector of size
#'          same as the number of locations in gwpca output or in a form of: 
#'          \code{sf$geomcolumn}.
#' @param method takes values either "membership" or "order". Is the method to 
#'               be used for grouping the spots together.
#' @gene.names a TRUE or FALSE value. Gives the option to translate the ENSGene 
#'             IDs to gene names. Defaults to FALSE. If TRUE, it MUST be 
#'             combined with a biomart table through the \code{biomart} argument.
#' @param biomart a biomart table that includes a "gene_name" column from ENSembl
#'                BioMart annotation tables. It can be generated using the 
#'                \code{biomaRt} package or the \code{create_biomart} function.
#' @param check.names a TRUE or FALSE value. It is passed to \code{data.frame}
#'                    inside \code{annotate_data}. Default is FALSE and if it is
#'                    TRUE it is only used when gene.names = TRUE as well. This
#'                    option is given because the check.names option can introduce
#'                    minor changes to the column names which then lead to errors
#'                    with the heatmap visualisation. A common change is "-" to 
#'                    ".". If check.names = FALSE then data.frame will not change
#'                    the column names. For more information look in 
#'                    \code{data.frame} help page.
#' 
#' 
#' @export

gwpca.topLeading.Gs <- function(gwpca, pc.no, genes.n, sf.geom, method = "membership",
                                gene.names = FALSE, biomart = NULL, 
                                check.names = FALSE){
    #### Check style of grouping ####
    if(method == "order"){
        message("Method selected for grouping is: order")
    } else if(method == "membership"){
        message("Method selected for grouping is: membership")
    } else {
        stop("Method provided is not valid.\nPlease select one of the methods: 'membership' OR 'order'")
    }
    
    #### Prepare data for the selected method ####
    message("Started grouping spots ...")
    local.loadings <- round(gwpca$loadings[,,pc.no], 4) # Retrieve the PC loadings
    dt <- apply(abs(local.loadings), 1, order, decreasing = TRUE) # Order the gene loadings on each spot in decreasing order
    dt_top <- dt[1:genes.n,]    # Select the top n genes with the highest loadings
    dt_top <- apply(dt_top, 2, function(x) colnames(local.loadings)[x]) # Add their names in the table
    if(gene.names){
        if(is.null(biomart)){
            stop("Please provide a biomart table")
        } else {
            message("Translating from ENSG IDs to Gene names ...")
            dt_top <- apply(dt_top, 2, function(x) 
                                annotate_vector(
                                    x,
                                    biomart,
                                    add = "gene_name",
                                    check.names = FALSE,
                                    output = "vector",
                                    annot.col = "gene_name"))
            message("Translation is done.")
            unique.n <- length(unique(c(dt_top))) # Count how many genes there are
            message("The number of individual genes found is: ", unique.n)
            message("The ", unique.n, " unique genes are:")
            print(table(c(dt_top)))
        }
    }
    
    #### Run the selected method ####
    if(method == "membership"){
        dt_top <- apply(dt_top, 2, sort) # Sort names alphabetically
        if(gene.names){ # If translation to gene names is active then this step is required before uniting the gene names
            dt_top <- bind_cols(dt_top)
        }
    }
    dt_top <- t(dt_top) %>%   # Concatenate gene names together for each spot
        as.data.frame() %>% 
        unite("Top_lead_Gs", 1:genes.n, sep = ";", remove = TRUE)
    group.n <- count(unique(dt_top)) # Count how many groups there are
    message("The number of individual groups found is: ", group.n) 
    
    #### Get them together in a df with geometries ready to plot ####
    dt_top <- data.frame(dt_top,
                         geometry = sf.geom) 
}
