#' @name add.perimeter.slide
#' @description This function will add a perimeter of off-tissue spots that
#' surround the Visium slide spots. This perimeter assists the Voronoi 
#' tessellation process by removing any instances where on-tissue spots at the 
#' edges end up without a polygon due to the way the tessellation is calculated.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @param data A data.frame as output from readSpacerangerMD.
#' @param spatialDir the pathway to the spatial directory that includes the 
#'                   .json file with the scale factors.
#'
#' @export

add.perimeter.slide <- function(data, spatialDir) {
    
    # Get the spot diameter and the max and min values for X and Y coordinates
    sp.diam <- spot_diameter(spatialDir)
    xmax <- round(max(data$pixel_x), 1)
    xmin <- round(min(data$pixel_x), 1)
    ymax <- round(max(data$pixel_y), 1)
    ymin <- round(min(data$pixel_y), 1)
    
    # Generate the perimeter spots
    dtXmax <- filter(data, pixel_x > (xmax - sp.diam)) %>% 
        mutate(pixel_x = pixel_x + 1.7*sp.diam)
    
    dtXmin <- filter(data, pixel_x < (xmin + sp.diam)) %>% 
        mutate(pixel_x = pixel_x - 1.7*sp.diam)
    
    dtYmax <- filter(data, pixel_y > (ymax - 2*sp.diam)) %>% 
        mutate(pixel_y = pixel_y + 2.8*sp.diam)
    
    dtYmin <- filter(data, pixel_y < (ymin + 2*sp.diam)) %>% 
        mutate(pixel_y = pixel_y - 2.8*sp.diam)
    
    dtPerim <- rbind(dtXmax, dtXmin, dtYmax, dtYmin) %>%
        mutate(Section = "perimeter") %>%
        mutate(Barcode = paste0(Barcode, ".perim"))
    
    return(rbind(data, dtPerim))
    
}

