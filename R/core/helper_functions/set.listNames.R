#' @name set.listNames
#' 
#' @description  A function to add names to list items from their attributes.
#' 
#' @param list the list to add names.
#' 

set.listNames <- function(list){
  n <- 1:length(list)
  
  named.list <- sapply(n, 
                       function(X){
                         names(list)[X] <- attr(list[[X]], "focus")
                         return(list[X])
                       })
  
  return(named.list)
}

