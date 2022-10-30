#' @name add.perimeter.tissue
#' @description This function will add a perimeter of off-tissue spots that are
#' neighbouring the on-tissue spots. All off-tissue spots are placed in bin 0
#' and all on-tissue spots in bin 1. The off-tissue perimeter spots are placed
#' in bin 2. This perimeter assists the Voronoi tessellation process by removing
#' any instances where non neighbouring spots would have been found as
#' neighbouring due to the way the tessellation is calculated.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @param data A data.frame as output from readSpaceranger.
#'
#' @export

add.perimeter.tissue <- function(data) {
    
    ## create a bin_0 subset
    bin_0 <- data %>%
        filter(Section == 0)
    
    ## create a bin_1 subset and add a new_bins column duplicate to Sections
    bin_1 <- data %>%
        filter(Section == 1) %>%
        mutate(new_bin = Section)
    
    ## find bin_0 spots with bin_1 neighbours and add them to bin_2.
    bin_0_2 <- cbind(bin_0, "new_bin" = apply(bin_0, 1, spot_neighbours, bin_1))
    
    return(rbind(bin_0_2, bin_1))

}
