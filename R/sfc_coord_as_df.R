#' @description This function converts an sfc object to a data.frame containing only the 
#' X and Y coordinates
#' 
#' 
#' @export


sfc_coord_as_df <- function(input) {
    ## convert to sf  get coordinates in a matrix and then convert to data.frame
    df <- input %>% 
        st_sf() %>% 
        st_coordinates() %>% 
        as.data.frame() %>%
        select(c("X", "Y"))
    
    ## Return
    df
}
