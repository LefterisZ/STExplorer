#' @name add.perimeter.slide
#' @description This function will add a perimeter of off-tissue spots that
#' surround the Visium slide spots. This perimeter assists the Voronoi 
#' tessellation process by removing any instances where on-tissue spots at the 
#' edges end up without a polygon due to the way the tessellation is calculated.
#' There have been cases where the pixel coordinates are a bit tilted towards 
#' one side and we need to adjust the perimeter addition. For that case we edit
#' the adjustselect and adjustadd parameters. 
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @param data A data.frame as output from readSpacerangerMD.
#' @param spatialDir the pathway to the spatial directory that includes the 
#'                   .json file with the scale factors.
#' @param adjustselect a vector of 4 numeric values to adjust the selection of
#'                     the on-the-edge slide spots for right, left, top and 
#'                     bottom (always in that order).
#' @param adjustadd a vector of 4 numeric values to adjust the addition of spots
#'                     on the perimeter for right, left, top and bottom (always 
#'                     in that order).
#'
#' @export

add.perimeter.slide <- function(data, spatialDir, 
                                adjustselect = c(1, 1, 2, 2),
                                adjustadd = c(1.7, 1.7, 2.8, 2.8)) {
    
    # Get the spot diameter and the max and min values for X and Y coordinates
    sp.diam <- spot_diameter(spatialDir)
    xmax <- round(max(data$pixel_x), 1)
    xmin <- round(min(data$pixel_x), 1)
    ymax <- round(max(data$pixel_y), 1)
    ymin <- round(min(data$pixel_y), 1)
    
    # Generate the perimeter spots
    dtXmax <- filter(data, pixel_x > (xmax - adjustselect[1]*sp.diam)) %>% 
        mutate(pixel_x = pixel_x + adjustadd[1]*sp.diam)
    
    dtXmin <- filter(data, pixel_x < (xmin + adjustselect[2]*sp.diam)) %>% 
        mutate(pixel_x = pixel_x - adjustadd[2]*sp.diam)
    
    dtYmax <- filter(data, pixel_y > (ymax - adjustselect[3]**sp.diam)) %>% 
        mutate(pixel_y = pixel_y + adjustadd[3]*sp.diam)
    
    dtYmin <- filter(data, pixel_y < (ymin + adjustselect[4]**sp.diam)) %>% 
        mutate(pixel_y = pixel_y - adjustadd[4]*sp.diam)
    
    dtPerim <- rbind(dtXmax, dtXmin, dtYmax, dtYmin) %>%
        mutate(Section = "perimeter") %>%
        mutate(Barcode = paste0(Barcode, ".perim"))
    
    return(rbind(data, dtPerim))
    
}

