#' This script contains test functions not implemented in STExplorer

test_apply <- function(x) {
    if (x["Section"] == 1) {
        return(1)
    } else {
        current_x <- as.integer(x["Spot_X"])
        current_y <- as.integer(x["Spot_Y"])
        neighbour_x <- c(current_x - 1, current_x + 1, 
                         current_x + 2, current_x - 2)
        neighbour_y <- c(current_y - 1, current_y, current_y + 1)
        neighbours <- test_tissue_positions[test_tissue_positions$Spot_X %in% neighbour_x & test_tissue_positions$Spot_Y %in% neighbour_y,]
        if (sum(neighbours$Section) > 0) {
            return(2)
        } else {
            return(0)
        }
    }
}
