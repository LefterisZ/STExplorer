#' @name get.nb.IDs
#' 
#' @description A function to recover neighbour IDs for a specific spot
#' 
#' 
#' @param obj a data.frame-like object that contains neighbours for each spot as
#'            a nested data frame.
#' @param barcode a character or a vector of characters with the Barcode names of 
#'             the spots we need to retrieve their neighbours.

get.nb.IDs <- function(obj, barcode){
    ## add a check here to check if there is a nb_IDs col available. If not look
    ## for a data col. If not again then print an error.
    for (b in 1:length(barcode)) {
        tmp_obj <- obj %>%
            filter(Barcode == barcode[b])
        
        nbs_df <- unnest(tmp_obj, cols = c("nb_IDs")) %>% # unnest the nested df
            select(c(nb_IDs, Barcode)) %>% # select only barcode columns
            st_drop_geometry() %>% # drop the geometry
            rename(is_nb_of = Barcode) %>% # rename the Barcode column
            left_join(obj, by = c("nb_IDs" = "Barcode")) %>% # pick nb info
            rbind(nbs_df) # rbind in case "barcode" vector has multiple barcodes
    }
    
    ## Return nbs
    return(nbs_df)
}