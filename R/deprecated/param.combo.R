#' @name param.combo
#' @description a GWPCA pipeline helper function to generate a parameters combo 
#'              df in order to run GWPCA with multiple combinations of parameters.
#' 
#' @param ... a series of variables and vectors in the same format as for a data 
#'            frame to be used for the combinations. The variable names will be
#'            used as column names.
#' @export

param.combo <- function(...){
    data <- expand.grid(..., stringsAsFactors = FALSE)
    data$kernl <- substring(data$kernel, 1, 3)
    data$obj <- paste0("pca_gw.", data$var, ".", data$k, ".", data$kernl)
    data$kernl <-NULL
    print(head(data, 3))
    
    return(data)
}