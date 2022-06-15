#' This script contains test functions not implemented in STExplorer

test_apply <- function(x, spots_all) {
    if (x["Section"] == 1) {
        return(1)
    } else {
        current_x <- as.integer(x["Spot_X"])
        current_y <- as.integer(x["Spot_Y"])
        neighbour_x <- c(current_x - 1, current_x + 1, 
                         current_x + 2, current_x - 2)
        neighbour_y <- c(current_y - 1, current_y, current_y + 1)
        neighbours <- spots_all[spots_all$Spot_X %in% neighbour_x & spots_all$Spot_Y %in% neighbour_y,]
        if (sum(neighbours$Section) > 0) {
            return(2)
        } else {
            return(0)
        }
    }
}

test_add_perimeter <- function(data) {
    
    tissue_spots_all <- data
    
    data_new <- cbind(data, "new_bin" = apply(data, 1, test_apply, tissue_spots_all))
    
    return(data_new)
}
