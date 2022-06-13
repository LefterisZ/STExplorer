#' @description This function takes the coordinates of the first spot of a 
#' clockwise sorted set of spots and adds it to the end. The input data is 
#' practically the output from \code{spot_sort()} function.
#' 
#' @param input an sfc object or the output from \code{spot_sort()}
#' @param return takes one of \code{"df"} or \code{"sfc"}. D
#' efault is \code{"df"}
#' 
#' 
#' 
#' @export

spot_add <- function(input, return = "df"){
    
    ## convert to data.frame
    df <- input %>%
        sfc_coord_as_df()
    
    ## get the first spot
    head_spot <- df %>%
        slice_head(n=1)
    
    ## add it at the bottom
    df <- df %>%
        add_row(head_spot, .after = nrow(df))

    ## Return
    if (return == "df") {
        return(df)
    } else if (return == "sfc") {
        return(st_as_sf(df, coords = c("X", "Y")) %>% 
                   st_combine())
    }
}
