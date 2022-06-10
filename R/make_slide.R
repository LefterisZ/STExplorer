#' @description This function generates a set of points that resemble the 10X Visium slides
#' The default will generate a slide of 5000 spots. You can adjust the slide
#' size to your liking.
#' 
#'Keep in mind that the Visium slide has an offset for both columns and rows. 
#'Meaning that each column X contains even only or odd only Y values. This 
#'generates the honeycomb-like structure of the 10X Visium slide
#' 
#' @param spots_x number of spots on the X axis (slide columns). 
#' Defaults to 100.
#' @param spots_y number of spots on the Y axis (slide rows). Defaults to 100.
#'  
#' 
#' @export

make_slide <- function(spots_x = 100, spots_y = 100) {
    
    ## Get slide X and Y ranges
    x <- 0:(spots_x - 1)
    y1 <- seq(0, (spots_y - 2), 2)
    y2 <- seq(1, (spots_y - 1), 2)
    
    ## combine x and y positions to generate a spot coordinates table
    df <- tibble(x = rep(x, each = spots_x/2),
                     y = rep(c(y1, y2), spots_y/2))
    
    
    ## Return
    return(df)
    
}
