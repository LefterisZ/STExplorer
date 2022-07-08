# This script will deal with extracting information needed for neighbour
#   identification and weighted matrix calculation from the Spaniel sce
#   object.

#packages needed
library(spdep)
library(sf)



# Export spot X-Y coordinates in a df
# xy_coord <- test_tissue_positions %>%
#   select(c("Barcode", "Spot_X", "Spot_Y", "Section")) %>%
#   remove_rownames() %>%
#   column_to_rownames(var = "Barcode") %>%
#   unite("X_Y", c(Spot_X, Spot_Y), remove = FALSE)

# function to get polygon from boundary box

bbox_polygon <- function(x) {
  bb <- sf::st_bbox(x)
  
  p <- matrix(
    c(bb["xmin"]-1, bb["ymin"]-1, 
      bb["xmin"]-1, bb["ymax"]+1,
      bb["xmax"]+1, bb["ymax"]+1, 
      bb["xmax"]+1, bb["ymin"]-1, 
      bb["xmin"]-1, bb["ymin"]-1),
    ncol = 2, byrow = T
  )
  
  sf::st_polygon(list(p))
}

# Prepare the dataset
mob_input <- seurat_list[["Olfactory_Bulb_A1_Results"]]@meta.data

mob_spot_position <- mob_input  %>% 
  select(c("Barcode", "Spot_X", "Spot_Y")) %>% 
  remove_rownames()

mob_centroids <- mob_spot_position %>% 
  st_as_sf(coords = c("Spot_X", "Spot_Y"))

# Each point has a separate geometry entry. 
head(mob_centroids)

# generate the hexagons
box <- st_sfc(bbox_polygon(mob_centroids))
head(box)

# This combines the points into a multipoint geometry:
mob <- st_union(mob_centroids)
head(mob)

# create a hexagon grid
grid <- st_make_grid(mob, 
                     cellsize = 1, 
                     square = FALSE,
                     offset = c(0, 0),
                     what = "polygons",
                     flat_topped = TRUE)
plot(grid, col = 0)

# Using the union of points generate a voronoi object
mob_voronoi <- st_voronoi(mob, box)
head(mob_voronoi)

# plot voronoi 
plot(mob_voronoi, col = 0) #rough
plot(st_intersection(st_cast(mob_voronoi), st_union(mob)), col = 1) # clip to smaller box
plot(mob_centroids, add = TRUE)


ggplot(xy_coord, aes(x = Spot_X, y = Spot_Y, label = X_Y, colour = Section)) +
  geom_point(size = 1) +
  xlim(15, 25) +
  ylim(30, 50) +
  geom_label()


ggsave("spatial_nb_coord_part.tiff",
       width = 834/96,
       height = 668/96,
       units = "in",
       dpi = 300)












plot(st_make_grid(mob_centroids, n = ), add = TRUE)




