#' @name readSpacerangerD (read Spaceranger Data)
#' @description This function takes the gene expression files 
#' from spaceranger and imports it as a data.frame for subsequent use.
#' 
#' @param dir The folder path where the gene expression count data are located.
#' 
#' 
#' @export

readSpacerangerD <- function(dir) {
    
    ## Read-in the gene counts as SingleCellExperiment (sce) object
    sce <- DropletUtils::read10xCounts(dir, version = "3")
    
    ## Get the gene counts as data.frame
    input_data <- data.frame(counts(sce))
    
    ## Add the barcodes of each spot as column names
    colnames(input_data) <- sce$Barcode
    
    ## Return
    return(input_data)
}

