#' @name spot_neighbours
#' @description This function will find which off-tissue spots (spots that have a
#' value of 0 in Section column) have an on-tissue spot neighbour. For every
#' spot that has a bin_1 neighbour returns the value of 2. For the rest returns
#' the value of 0.
#'
#' This is a function that will work only inside an \code{apply()} function
#' combined with a cbind (see example).
#'
#' @param row This is the x argument of the \code{FUN(X)} that \code{apply()}
#'           will apply on each row of the \code{apply()} function's input.
#'
#' @param bin.data This is a data.frame in which the function will try to identify
#'                if there are neighbours of the spot in question. This
#'                data.frame contains the positions of all the on-tissue spots
#'                (bin 1)
#'
#' @export
#' @example bin_0_and_2 <- cbind(bin_0, "new_bin" = apply(bin_0, 1, spot_neighbours, bin_1))

spot_neighbours <- function(row, bin.data) {
    curr_X <- as.integer(row["Spot_X"]) # get current X
    curr_Y <- as.integer(row["Spot_Y"]) # get current Y
    
    ## calculate the putative neighbour (nb) X and Y positions.
    nb_y <- c(curr_Y - 1, curr_Y, curr_Y + 1)
    nb_x <- c(curr_X - 1, curr_X + 1,
              curr_X + 2, curr_X - 2)
    
    ## find neighbours (nbs) from nb inside bin_1
    nbs <-
        bin.data[bin.data[["Spot_X"]] %in% nb_x &
                     bin.data[["Spot_Y"]] %in% nb_y,]
    
    ## if bin_0 spot has 1 or more neighbours in bin_1 return 2 (for bin_2)
    ## otherwise return 0
    if (sum(nbs[["Section"]]) > 0) {
        return(2)
    } else {
        return(0)
    }
    
}
