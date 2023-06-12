#' @name add.spotHex
#' @description This function will add Hexagons in the SFE object utilising 
#' a perimeter of off-tissue spots that
#' surround the Visium slide spots. This perimeter assists the Voronoi 
#' tessellation process by removing any instances where on-tissue spots at the 
#' edges end up without a polygon due to the way the tessellation is calculated.
#' There have been cases where the pixel coordinates are a bit tilted towards 
#' one side and we need to adjust the perimeter addition. For that case we edit
#' the adjustselect and adjustadd parameters. 
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @param obj The SpatialFeaturesExperiment object.
#' @param samples 	a character vector specifying one or more directories, each 
#' corresponding to a 10x Genomics Visium sample (see Details); if provided, 
#' names will be used as sample identifiers.
#' @param sample_id character string specifying unique sample identifiers, one 
#' for each directory specified via samples; ignored if 
#' !is.null(names(samples)).
#' @param res the desired resolution. Can take one of "lowres", "hires", 
#' "fullres".
#' @export

# The idea is to work on the array rows and columns. If there are spots on the 
# first and last row/ column then we select those, make a TRUE/FALSE vector to 
# fetch the coordinates from the spatialCoordinates. Then add the diameter*1.5, 
# then take their new coordinates and add them back to spatialCoordinates 
# (problems with colData and SFE? Maybe export spatialCoordiantes work and then 
# add to SFE ONLY the hexagons.)

add.spotHex <- function(obj,
                        samples,
                        sample_id,
                        res = c("lowres", "hires", "fullres")) {
    ## Prepare required data
    res <- match.arg(res)
    res <- .int_resSwitch(res)
    cData <- cbind(colData(obj), spatialCoords(obj))
    n <- length(samples)
    dataList <- vector("list", length = n)
    
    for (i in seq_along(samples)) {
        ## Get spot diameter
        .sp_diam <- metadata(obj)$spotDiameter[[sample_id[i]]][[res]]
        data <- read.csv(file.path(samples[i], "outs/spatial", 
                                   "tissue_positions_list.csv"), 
                         stringsAsFactors = FALSE, header = FALSE)
        colnames(data) <- c("Barcode", "Section", "Spot_Y", 
                            "Spot_X", "Image_Y", "Image_X")
        ## Get min/max values for the first/last capture area rows/columns. It 
        ## is not fixed because Visium can have 6.5mm or 11mm capture areas.
        int_list_minMax <- .int_spotHex_minMax(data)
        
        
        ## Fetch data and coordinates for the required locations
        int_list_subset <- .int_spotHex_subset(int_list_minMax, cData, data)
        
        ## Generate the perimeter spots
        int_df_perim <- .int_spotHex_gen(int_list_subset, .sp_diam)
        
        ## Add them to the rest of the data
        dtPerim <- rbind(int_df_perim,
                         data)
        
        ## Convert spots to centroids
        centroids <- as.data.frame(dtPerim) %>% 
            st_as_sf(coords = c("Image_X", "Image_Y"), 
                     remove = TRUE)
        
        ## Combine the points into a multipoint geometry:
        cntd_union <- st_union(centroids)
        
        ## Use the union of points to generate a voronoi object
        voronoi <- st_voronoi(cntd_union, bOnlyEdges = TRUE)
        
        ## Create an enveloped voronoi tessellation around the tissue
        voronoi_env <- st_intersection(st_cast(voronoi), 
                                       st_convex_hull(cntd_union))
        
        ## Generate the POLYGONS from the MULTILINESTRING
        polygons <- st_polygonize(voronoi_env) %>% # polygonise the tessellation
            st_cast() %>% # convert GEOMETRYCOLLECTION to multiple POLYGONS
            st_sf() %>%  # convert sfc object to sf for st_join afterwards
            st_join(., 
                    centroids[centroids$Section == 1,],
                    join = st_contains,
                    left = FALSE) %>% # Join the centroids with the POLYGONS
            arrange(Barcode) %>% 
            dplyr::select(geometry)
        
        ## Append it to a list
        dataList[[i]] <- polygons
    }
    
    hexes <- dplyr::bind_rows(dataList)
    rownames(hexes) <- colnames(obj)
    colGeometry(obj, "spotHex") <- hexes
    
    return(obj)
    
}

