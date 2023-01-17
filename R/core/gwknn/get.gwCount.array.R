#' @name get.gwCount.array
#' 
#' @description A function to produce a 3D array with dimensions s x g x s. 
#'              Where s = spots (locations) and g = genes. Takes as input an 
#'              s x g matrix of gene counts, an s x s distance matrix with 
#'              weighted distances between locations. Returns as output a 3D
#'              array with gene expression counts weighted for each location (or
#'              for a number of selected locations).
#' 
#' @param obs A matrix containing the observation data, usually gene expression 
#'            data where columns are the genes and rows are the locations.
#' @param wdmat An n x n matrix (n = number of locations) of weighted distances 
#'              between all locations.
#' @param focus.n Numeric or character. The indexes (numeric) or the names 
#'                (character) of the locations you want in focus. It can be a 
#'                vector of indexes from locations that are of interest. Default
#'                behaviour is for focus.n to be missing. This will result to  
#'                all locations being considered.
#' 
#' @return A 3D array with dims = [s, g, focus.n]
#' 
#' @export


get.gwCount.array <- function(obs, wdmat, focus.n){
    
    # Check indexes of locations 
    if(missing(focus.n)){
        focus.n <- 1:nrow(obs) # indexes of locations to use (z-axis)
    } else if(is.vector(focus.n) & is.numeric(focus.n)){
        message("A selection of location indexes was provided...")
        message("Locations with indexes: ", paste(focus.n, collpse = " "))
    } else if(is.vector(focus.n) & is.character(focus.n)) {
        message("A selection of location names was provided...")
        message("Locations with names: ", paste(focus.n, collpse = " "))
        focus.n <- match(focus.n, dimnames(obs)[[1]])
    } else {
        stop("The vector provided at the focus.n argument does not contain
             numeric values OR character values only. 
             Please make sure you provide location indexes only.")
    }
    
    # Weight expression data
    temp <- sapply(focus.n, .get.gw.counts, 
                   obs = obs, wdmat = wdmat, 
                   simplify = "array")
    
    ## add dimnames --> spot names as rows and columns
    dimnames(temp)[[1]] <- rownames(obs)
    dimnames(temp)[[2]] <- colnames(obs)
    dimnames(temp)[[3]] <- rownames(obs)[focus.n]
    
    return(temp)
    
}
