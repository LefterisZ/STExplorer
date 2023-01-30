#' @name .dist.wrap
#' 
#' @description a wrapper for the base dist function to be used in a future_sapply
#'              function
#' 
#' 
#' 

.dist.wrap <- function(X, method, p, ...){
    d <- dist(gwCounts[,,X], method, p) %>% as.matrix()
    Sys.sleep(0.001)
    pr()
    
    return(d)
}

