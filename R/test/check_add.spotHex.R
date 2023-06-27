################################################################################
# This script checks the outputed hexagons.
# Use it when the `addGeometries()` function returns error for Hexagons.
# First add the `assign()` functions inside the loop of `add.spotHex()` apart 
#   from the `assign(hexes,...)` which should be adde dfter the `hexes` object 
#   is created outside of the loop.
# Use the `ggplot()` to visualise.
# Use the `rm()` to clear the environment.

assign("cData", cData, envir = .GlobalEnv)
assign("sp_diam", .sp_diam, envir = .GlobalEnv)
assign("data", data, envir = .GlobalEnv)
assign("int_list_minMax", int_list_minMax, envir = .GlobalEnv)
assign("int_list_subset", int_list_subset, envir = .GlobalEnv)
assign("int_df_perim", int_df_perim, envir = .GlobalEnv)
assign("dtPerim", dtPerim, envir = .GlobalEnv)
assign("centroids", centroids, envir = .GlobalEnv)
assign("cntd_union", cntd_union, envir = .GlobalEnv)
assign("voronoi", voronoi, envir = .GlobalEnv)
assign("voronoi_env", voronoi_env, envir = .GlobalEnv)
assign("polygons", polygons, envir = .GlobalEnv)
assign("dataList", hexes, envir = .GlobalEnv)

assign("hexes", hexes, envir = .GlobalEnv)

# ---------------------------------------------------------------------------- #
array_sizeX <- data %>%
    filter(Spot_X == 0) %>% 
    select(Image_Y, Image_X)
array_sizeY <- data %>%
    filter(Spot_Y == 0) %>% 
    select(Image_Y, Image_X)

ggplot() + 
    geom_sf(data = voronoi_env) + 
    geom_sf(data = centroids, aes(geometry = geometry)) + 
    geom_sf(data = hexes, aes(geometry = geometry)) + 
    geom_point(data = as.data.frame(spatialCoords(sfe)), aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres), colour = "red") + 
    geom_point(data = array_sizeX, aes(x = Image_X, y = Image_Y), colour = "green") + 
    geom_point(data = array_sizeY, aes(x = Image_X, y = Image_Y), colour = "blue")


rm(hexes,sp_diam,data,int_list_minMax,int_list_subset,int_df_perim,dtPerim,
   centroids,cntd_union,voronoi,voronoi_env,polygons,dataList,array_sizeX, 
   array_sizeY,cData)

