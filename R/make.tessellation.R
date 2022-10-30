#' @name  make.tessellation
#' 
#' @description 
#' 
#' @param 
#' 
#' @export

make.tessellation <- function(input){
    ## Select spots in both bins (Sections) 0 and 1 ----
    spot_position <- input %>% 
        select(c("Barcode", "pixel_x", "pixel_y", "Section"))
    
    ## Convert spots to centroids ----
    centroids <- spot_position %>% 
        st_as_sf(coords = c("pixel_x", "pixel_y"), 
                 remove = FALSE)
    
    ## Combine the points into a multipoint geometry: ----
    cntd_union <- st_union(centroids)
    head(cntd_union)
    
    ## Use the union of points to generate a voronoi object ----
    voronoi <- st_voronoi(cntd_union, bOnlyEdges = TRUE)
    head(voronoi)
    
    ## Create an enveloped voronoi tessellation around the tissue ----
    voronoi_env <- st_intersection(st_cast(voronoi), st_convex_hull(cntd_union))
    head(voronoi_env)
    
    return(voronoi_env)
    
}