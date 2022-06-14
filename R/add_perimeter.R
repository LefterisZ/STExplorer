#' @description This function will add a perimeter of off-tissue spots that are
#' neighbouring the on-tissue spots. All off-tissue spots are placed in bin 0
#' and all on-tissue spots in bin 1. The off-tissue perimeter spots are placed
#' in bin 2. This perimeter assists the Voronoi tessellation process by removing
#' any instances where non neighbouring spots would have been found as
#' neighbouring due to the way the tessellation is calculated.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @param positions A data.frame as output from readSpaceranger.
#' @param
#'
#' @export

add_perimeter <- function(data) {
    ## add a bin column to contain the new bins and populate it with the already
    ## existing ones from "Section" column.
    positions <- data %>%
        mutate(bins = Section)
    
    ## create a bin_0 subset
    bin_0 <- positions %>%
        filter(Section == 0)
    
    bin_1 <- positions %>%
        filter(Section == 1)
    
    ## find bin_0 spots with bin_1 neighbours and add them to bin_2.
    x <- apply(bin_0, 1, spot_neighbours, bin.1 = bin_1)
    
    return(as.data.frame(rbind(bin_0, bin_1)))

}
