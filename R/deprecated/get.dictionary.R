#' @name get.dictionary
#' 
#' @description A function to create a dictionary of matches between the 
#'              spot Barcodes and the rownames as well as the gene names and the
#'              column names.
#' 
#' @param obs the gene expression matrix with Barcodes as row names and ENSGIDs 
#'            as column names
#' 
#' @return a list of row name and column name matches to numerical row and 
#'         column names.
#' 
#' @export

get.dictionary <- function(obs){
  ## Sort rownames and colnames
  obs <- obs %>%
    .[order(rownames(.)),] %>% 
    .[,order(colnames(.))]
  
  ## Create spots dictionary
  spots_dict <- as.matrix(rownames(obs))
  rownames(spots_dict) <- 1:nrow(obs)
  colnames(spots_dict) <- "Barcode"
  
  ## Create genes dictionary
  genes_dict <- as.matrix(colnames(obs))
  rownames(genes_dict) <- paste0("g",1:ncol(obs))
  colnames(genes_dict) <- "ENSG_ID"
  
  ## Combine both to a list
  dict <- list("spots" = spots_dict,
               "genes" = genes_dict)
  
  return(dict)
  
}
