#' @name gwpca.Leading.G.single 
#' 
#' @param gwpca a gwpca output.
#' @param pc.no the Principal Component (PC) to be examined.
#' @param genes.n an integer indicating how many genes you want to be included.
#' @sf.geom the geometry column from an sf object containing a polygons 
#'          geometry column to be used for mapping. It can be a vector of size
#'          same as the number of locations in gwpca output or in a form of: 
#'          \code{sf$geomcolumn}.
#' @gene.names a TRUE or FALSE value. Gives the option to translate the ENSGene 
#'             IDs to gene names. Defaults to FALSE. If TRUE, it MUST be 
#'             combined with a biomart table through the \code{biomart} argument.
#' @param biomart a biomart table that includes a "gene_name" column from ENSembl
#'                BioMart annotation tables. It can be generated using the 
#'                \code{biomaRt} package or the \code{create_biomart} function.
#' @param check.names a TRUE or FALSE value. It is passed to \code{data.frame}
#'                    inside \code{annotate_data.frame}. Default is FALSE and if 
#'                    it is TRUE it is only used when gene.names = TRUE as well. 
#'                    This option is given because the check.names option can 
#'                    introduce minor changes to the column names which then 
#'                    lead to errors with the heatmap visualisation. A common 
#'                    change is "-" to ".". If check.names = FALSE then 
#'                    data.frame will not change the column names. For more 
#'                    information look in \code{data.frame} help page.
#' 
#' 
#' @export

gwpca.Leading.G.single <- function(gwpca, pc.no, sf.geom,
                                   gene.names = FALSE, biomart = NULL,
                                   check.names = FALSE){
    # Create the PC name
    pc.name <- paste0("PC", pc.no)
    
    # Get local loadings for pc.no
    local.loadings <- round(gwpca$loadings[,,pc.no], 4) 
    
    # Find the leading item (colname) on each location (rows)
    message("Finding leading genes ...")
    lead.item <- colnames(local.loadings)[max.col(abs(local.loadings))] %>%
        as.data.frame()
    colnames(lead.item) <- pc.name
    
    # Translate to gene names
    if(gene.names){
        if(is.null(biomart)){
            stop("Please provide a biomart table")
        } else {
            message("Translating from ENSG IDs to Gene names ...")
            lead.item <- annotate_data.frame(
                    data = lead.item,
                    biomart,
                    add = "gene_name",
                    check.names = FALSE,
                    column = 1,
                    output = "data.frame")
            message("Translation is done.")
            unique.n <- length(unique(c(lead.item$gene_name))) # Count how many genes there are
            message("The number of leading genes found is: ", unique.n)
        }
    }
    
    # Print how many items are leading in this PC
    message("The leading genes in ", pc.name, " are: ") 
    print(table(lead.item$gene_name))
    
    # Make the outputted data frame
    lead.item <- data.frame(lead.item,
                            sf.geom) %>% 
        dplyr::select( -all_of(c("id")))
    
    colnames(lead.item) <- c(pc.name, "ENSG_ID", "geometry")
    
    return(lead.item)
}
