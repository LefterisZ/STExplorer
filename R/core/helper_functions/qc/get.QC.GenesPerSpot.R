#' @name get.QC.GenesPerSpot
#' 
#' @description get the number of genes found in each spot
#' 
#' @param count_table a gene counts table with genes in rows and spots in columns
#' 
#' @param select defaults to NULL and proceeds to work on the complete table. 
#'               A vector of columns can be given to select them and work on this
#'               selection only. 
#' 
#' @export

get.QC.GenesPerSpot <- function(count_table, select = NULL) {
  
  ## Select specific spots
  if (!is.null(select)) {
    tmp <- count_table %>% 
      select(all_of(select))
  } else {
    tmp <- count_table
  }
  
  ## Find how many genes are expressed (gene expression > 0)
  tmp <- tmp > 0
  tmp <- as.data.frame(tmp) %>% 
    colSums() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Barcode") %>%
    rename("gene_number" = ".")
  
  return(tmp)
  
}
