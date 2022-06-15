#'@description This function will find the neighbours of each spot
#'
#'
#'
#'@export

spot_neighbours <- function(bin.0_row, bin.1) {
    
    curr_X <- as.integer(bin.0_row["Spot_X"]) # get current X
    curr_Y <- as.integer(bin.0_row["Spot_Y"]) # get current Y
    
    ## calculate the putative neighbour (nb) X and Y positions.
    nb_x <- c(curr_X - 1, curr_X, curr_X + 1)
    nb_y <- c(curr_Y - 1, curr_Y + 1, 
              curr_Y + 2, curr_Y - 2)
    
    ## find neighbours (nbs) from nb inside bin_1
    nbs <- bin.1[bin.1$Spot_X %in% nb_x & bin.1$Spot_Y %in% nb_y,]
    
    ## if bin_0 spot has 1 or more neighbours in bin_1 return 2 (for bin_2)
    ## otherwise return 0
    if (sum(nbs[["Section"]]) >0 ){
        return(2)
    } else {
        return(0)
    }
    
}
