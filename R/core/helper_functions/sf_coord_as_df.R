#' @description This function takes an sf object and extracts the X and Y coordinates as a 
#' data.frame
#' 
#' 
#' @export

sf_coord_as_df <- function(input_sf) {
    
    ## convert to data.frame (make it a separate function)
    df <- input_sf %>% 
        st_coordinates() %>% 
        as.data.frame(input_sf) %>%
        select(c("X", "Y"))
    
    ## Return
    return(df)
}
