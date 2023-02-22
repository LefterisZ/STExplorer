#' a function to update an existing object
#' 

`updt<-` <- function(x, ..., value){
  ## x is the object to be manipulated, value the object to be assigned
  x[... ,] <- value
  x
}
