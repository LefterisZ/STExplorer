#' @name spot_sort
#' 
#' @description  Function to sort a set of spots in a clockwise or 
#' anti-clockwise direction. These spots can be randomly selected spots from 
#' the 10X Visium slide layout. Each spot has a set of coordinates that this 
#' function will sort.
#' 
#' The function calculates the position of a new central spot using the median 
#' for each spot's coordinates (median of all X positions, median of all Y 
#' positions). Then, uses this central spot coordinates to calculate the angles 
#' for the other spots. These angles are then used to sort the spots in a
#' clockwise/ or anti-clockwise way. 
#' 
#' \code{sort_spots} is useful to generate a random perimeter that resembles
#' a tissue section. This is also useful to simulate random spatial 
#' transcriptomics datasets
#' 
#' @author Eleftherios (Lefteris) Zormpas
#' 
#' @param input sf or sfc class object containing a MULTIPOINT set of x and y 
#' coordinates for the selected spots. 
#' @param x Name of x variable in \code{input}. Defaults to "X".
#' @param y Name of y variable in \code{input}. Defaults to "Y". 
#' @param clockwise If \code{TRUE}, order spots clockwise, If \code{FALSE}, 
#' order them anti-clockwise. 
#'
#' @returns an sfc object
#'
#' @export


spot_sort <- function(input, x = "X", y = "Y", clockwise = TRUE) {
    
    ## check class if sf, or sfc or data.frame, if sf or sfc make them data.frame
    if (is(input, "sf")) {
        
        ## convert to data.frame
        df <- sf_coord_as_df(input)
        
    } else if (is(input, "sfc")) {
        
        ## convert to sf  get the coordinates and then convert to data.frame
        df <- sfc_coord_as_df(input)
        
    } else if (is(input, "data.frame")) {
        
        ## notify that input is already of class data.frame
        message("The input is already of class data.frame. No conversion is needed.")
    }
    
    ## Calculate the central point of the imported points
    ##(make it a separate function)
    x_central <- median(df[, x])
    y_central <- median(df[, y])
    
    ## Calculate the angles and sort
    df <- df %>%
        ## calculate differences of each spot from the central one
        mutate(diff_x = X - x_central, diff_y = Y - y_central) %>% 
        ## calculate the angles
        mutate(angles = atan2(diff_y, diff_x)) %>% 
        ## sort rows
        {if (clockwise) {
            arrange(., desc(angles))
        } else {arrange(., angles)}} %>% 
        ## remove redundant columns
        select(c(c("X", "Y")))
    
    
    ## Return
    return(st_as_sf(df, coords = c("X", "Y")) %>%
               st_combine())
    
}
