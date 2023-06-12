#' @name get.QC.CountsPerSpot
#' 
#' @description get the total expression of all genes per spot
#' 
#' @param count_table a gene counts table with genes in rows and spots in columns
#' 
#' @param select defaults to NULL and proceeds to work on the complete table. 
#'               A vector of columns can be given to select them and work on this
#'               selection only. 
#' 
#' @export

get.QC.CountsPerSpot <- function(count_table, select = NULL) {
  
  ## Select specific spots
  if (!is.null(select)) {
    tmp <- count_table %>% 
      select(all_of(select))
  } else {
    tmp <- count_table
  }
  
  ## Get the total gene expression per spot
  tmp <- tmp %>% 
      colSums() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "Barcode") %>%
      dplyr::rename("CountsPerSpot" = ".") %>% 
      arrange(CountsPerSpot) %>%
      mutate(aa = 1:nrow(.),
             log2_CountsPerSpot = log2(CountsPerSpot),
             log1p_CountsPerSpot = log1p(CountsPerSpot))
  
  knee <- kneedle(tmp$aa, tmp$CountsPerSpot)
  attr(tmp, "knee") <- knee
  
  return(tmp)
  
}
