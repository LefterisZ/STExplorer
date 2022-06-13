#'@description This function will find the neighbours of each spot
#'
#'
#'
#'@export

spot_neighbours <- function(input) {
    x <- input["Spot_X"]
    y <- input["Spot_Y"]
    neighbours <- cbind(nb_x = c(x - 1, x, x + 1, x + 1, x, x - 1),
                        nv_y = c(y - 1, y - 2, y - 1, y + 1, y + 2, y + 1))
}
