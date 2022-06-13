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

add_perimenter <- function(data) {
    ## add a bin column to contain the new bins and populate it with the already
    ## existing ones from "Section" column.
    positions <- data %>% 
        mutate(bins = Section)
    
    ## create a bin_0 subset
    bin_0 <- positions %>%
        filter(Section == 0)
    
    ## find neighbouring spots and their bin. add the bin 0 spots into bin 2.
    group_modify(bin_0, )
    
    ## pick spots from bin 2 and add them into the tessellation input
    
    
}
