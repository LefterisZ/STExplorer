#'@description This function will find the neighbours of each spot
#'
#'
#'
#'@export

spot_neighbours <- function(bin.0_row, bin.0, bin.1) {

    Spot_X <- as.integer(bin.0_row["Spot_X"]) #get current X
    Spot_Y <- as.integer(bin.0_row["Spot_Y"]) #get current Y
    
    ## c-bind them and convert them to a data.table. Add putative bins = 2.
    cur_XY <- cbind(Spot_X = Spot_X, Spot_Y = Spot_Y, bins = 2) %>%
        as.data.frame() %>% 
        data.table::setDT()
    
    ## calculate the 6 putative neighbours (nb) and convert to data.table.
    nb <- cbind(Spot_X = c(Spot_X - 1, Spot_X, Spot_X + 1, 
                           Spot_X + 1, Spot_X, Spot_X - 1), 
                Spot_Y = c(Spot_Y - 1, Spot_Y - 2, Spot_Y - 1, 
                           Spot_Y + 1, Spot_Y + 2, Spot_Y + 1), 
                bins = rep(2, 6)) %>% 
        as.data.frame() %>% 
        data.table::setDT()
    
    ## find neighbours from nb inside bin_1 and find if bin sum is > 0. If yes
    ## change bin number in bin_0 to 2.
    if (sum(bin_1[nb, on = c("Spot_X", "Spot_Y"), bins], na.rm = TRUE) >0 ){
        bin.0_row["bins"] = 2
    } 
    
}
