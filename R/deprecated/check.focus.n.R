#' @name check.focus.n
#' 
#' @description a function to check the input vector of locations in focus
#' 
#' @param .focus.n a character vector with the indexes from the location(s) in 
#'                focus.
#' @param .dict the dictionary that matches indexes with barcode names.
#' 
#' 

check.focus.n <- function(.focus.n = NULL, .dict){
  
  # Check indexes of locations
  if (is.null(.focus.n)) {
    .focus.n <- rownames(.dict[["spots"]]) # indexes of locations to use
    
  } else if (is.vector(.focus.n)) {
    vec <- as.numeric(.focus.n)
    m <- match(.focus.n, rownames(.dict[["spots"]]))
    sum.vec <- sum(is.na(vec))
    sum.m <- sum(is.na(m))
    l <- length(.focus.n)
    
    ## Check if the vector is with numbers as characters only or barcodes only
    if (sum.vec > 0 & sum.vec < l) {
      stop("The vector provided at the focus.n argument does not contain
            character values only. 
            Please make sure you provide EITHER location indexes as characters 
            OR barcode names only. NEITHER integers NOR mixture of indexes and
            barcodes")
      
    } else if (sum.m != 0) {
      message("A selection of location names was provided...")
      message("Locations with names: ", paste(.focus.n, collpse = " "))
      .focus.n <- match(.focus.n, .dict[["spots"]]) # find the indexes
    
    } else {
      message("A selection of location indexes was provided...")
      message("Locations with indexes: ", paste(.focus.n, collpse = " "))
      .focus.n <- .focus.n
    }
  } else {
    stop("Please provide a valid argument for the locations in focus. Either 
         leave it as NULL or provide a character vector.")
  }
  
  return(.focus.n)
}
