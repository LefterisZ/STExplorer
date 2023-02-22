#' @name data.normalise
#' 
#' @description A function to normalise the gene expression counts table
#' 
#' @param data the counts table
#' 
#' @export
#' 

data.normalise <- function(data){

    tmp <- data + 1
    
    return(
      log2(x = scale(x = tmp, center = FALSE, scale = colSums(x = tmp)) * 10000)
    )
  
}
