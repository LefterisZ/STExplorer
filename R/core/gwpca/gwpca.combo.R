#' @name gwpca.combo
#' @description a GWPCA pipeline wrapper function to run multile gwpcas and  
#'              and output the results to a list. The function engulfs all the 
#'              rest of the multiple gwpca helper functions: \code{combo.wrap()}
#'              that uses data generated from \code{param.combo()} to pass them
#'              to \code{gwpca.param.combo()}.
#' @param dt.combo the df with the parameter combinations
#' 
#' @export

# function that wraps gwpca.param.combo and outputs objects in a list
gwpca.combo <- function(dt.combo){
    
    pcagw.list <- list()
    
    pcagw.list <- with(dt.combo,
                        mapply(combo.wrap, 
                               var.no, k, kernel, 
                               MoreArgs = list(list = pcagw.list)))
    
    pcagw.list <- setNames(pcagw.list, dt.combo$obj)
    
    return(pcagw.list)
}