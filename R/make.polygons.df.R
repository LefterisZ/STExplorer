#' @name make.polygons.df
#' 
#' @description 
#' 
#' @param 
#' 
#' 
#' @export

#' @note this doesn't include make neighbours. Neighbours list needs to be added
#'       as input.

make.polygons.df <- function(){
    ## Generate the POLYGONS from the MULTILINESTRING for bin_1 only and attach the  
    ##  barcode names
    polygons <- st_polygonize(voronoi_env) %>% # polygonise the tessellation
        st_cast() %>% # convert GEOMETRYCOLLECTION to multiple POLYGONS
        st_sf() %>%  # convert sfc object to sf for st_join afterwards
        st_join(., 
                centroids[centroids$Section == 1,],
                join = st_contains,
                left = FALSE) %>% # Join the centroids with the POLYGONS
        mutate(Barcode_rn = Barcode) %>% # duplicate the barcode column
        column_to_rownames("Barcode_rn") %>% # move duplicate column to row names
        st_sf() # convert back to sf (mutate makes it a df)
    
    ## Add number of neighbours for each polygon back to the polygons object
    polygons$nb_count <- card(neighbours)
    
    ## Add the neighbour (nb) IDs as a nested df in the polygons object
    nb_IDs <- neighbours %>%
        nb2lines(., coords = polygons$geometry) %>% #get nb connecting lines 
        as("sf") %>% #convert to sf
        st_drop_geometry() %>% #drop geometry column
        select(i_ID, j_ID) %>% #select only nb ID columns
        rename(nb_IDs = j_ID) %>% #rename the neighbours ID column
        group_by(i_ID) %>% #group by spot
        nest() #nest the groupings
    
    polygons <- right_join(polygons, nb_IDs, by = c("Barcode" = "i_ID")) %>%
        rename(nb_IDs = data, geom_pol = geometry)
    
    
    ## Update the polygon object to keep the centroid geometries as well
    polygons <- left_join(as.data.frame(polygons), as.data.frame(centroids), 
                          by = c("Barcode" = "Barcode"), suffix = c("", ".y")) %>%
        select(!ends_with(".y")) %>% 
        rename(geom_cntd = geometry) %>% 
        st_sf(sf_column_name = "geom_pol")
    
}