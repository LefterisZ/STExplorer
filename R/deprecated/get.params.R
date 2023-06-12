#' @name get.params
#' @description a GWPCA pipeline helper function to extract the parameters from
#'              every gwpca output object. The function is ideal to be used with
#'              sapply over a list of multiple gwpca objects.
#' @param gwpca.out a gwpca output object
#' 
#' @export

get.params <- function(gwpca.out){
    obj <- data.frame("vars" = length(gwpca.out$GW.arguments$vars),
                      "spots" = gwpca.out$GW.arguments$dp.n,
                      "PCs" = gwpca.out$GW.arguments$k,
                      "kernel" = gwpca.out$GW.arguments$kernel,
                      "minutes" = difftime(gwpca.out$timings$stop, 
                                           gwpca.out$timings$start,
                                           units = "mins"))
    
    return(obj)
}
