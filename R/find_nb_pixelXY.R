# This script will deal with extracting information needed for neighbour
#   identification and weighted matrix calculation from the Spaniel sce
#   object.

#packages needed
library(spdep)
library(sf)
library(jsonlite)
library(ggplot2)

source(file = file.path(scriptsDir, "spot_diameter.R"))
source(file = file.path(scriptsDir, "make_bb_polygon.R"))




# Export spot X-Y coordinates in a df
# xy_coord <- seurat_list[["Olfactory_Bulb_A1_Results"]]@meta.data %>%
#   select(c("Barcode", "Spot_X", "Spot_Y")) %>% 
#   remove_rownames() %>% 
#   column_to_rownames(var = "Barcode") %>% 
#   unite("X_Y", c(Spot_X, Spot_Y), remove = FALSE)

# Prepare the dataset
mob_input <- seurat_list[["Olfactory_Bulb_A1_Results"]]@meta.data

mob_spot_position <- mob_input  %>% 
  select(c("Barcode", "pixel_x", "pixel_y")) %>% 
  remove_rownames() 

mob_centroids <- mob_spot_position %>% 
  st_as_sf(coords = c("pixel_x", "pixel_y"))

# Each point has a separate geometry entry. 
head(mob_centroids)

# generate the bounding box

  # get the folder path to the scale factors .json file
spatialDir <- file.path(dataDir, sampleInfo$fileFolders, "spatial")

  # calculate spot diameter
spot_diam <- spot_diameter(spatialDir, "scalefactors_json.json")

  # Get a polygon from boundary box
box <- st_sfc(make_bb_polygon(mob_centroids, spot_diam))
head(box)

boxXmin <- min(box[[1]][[1]][,1])
boxXmax <- max(box[[1]][[1]][,1])
boxYmin <- min(box[[1]][[1]][,2])
boxYmax <- max(box[[1]][[1]][,2])


# This combines the points into a multipoint geometry:
mob <- st_union(mob_centroids)
head(mob)

# Using the union of points generate a voronoi object
mob_voronoi <- st_voronoi(mob, bOnlyEdges = FALSE)
head(mob_voronoi)

# get the line coordinates, filter them for X-Y min and X-Y max and plot them?
# intersect the mob_voronoi with the convex hull.
# how to find the surrounding neighbours

# plot voronoi 
plot(mob_voronoi, col = 0, axes = TRUE, xlim = c(boxXmin, boxXmax), ylim = c(boxYmin, boxYmax)) #rough
plot(mob_centroids, col = "red", add = TRUE)
plot(st_convex_hull(mob), border = "blue", col = rgb(1, 1, 1, 0.0), add = TRUE)

# plot voronoi enveloped
mob_voronoi_env <- st_intersection(st_cast(mob_voronoi), st_convex_hull(mob))
plot(mob_voronoi_env, col = 0, axes = TRUE)
plot(mob_centroids, col = "red", add = TRUE)
plot(st_convex_hull(mob), border = "blue", col = rgb(1, 1, 1, 0.0), add = TRUE)

#plot test voronoi
plot(test_mob_voronoi, col = 0, axes = TRUE) #rough


ggplot() +
  geom_sf(data = mob_voronoi_env, colour = "black", fill = "white") + 
  geom_sf(data = mob_centroids, colour = "red3") + 
  xlim(boxXmin, boxXmax) + 
  ylim(boxYmin, boxYmax) + 
  # Add titles and visually format the plot:
  labs(title = paste("Voronoi tessellation"),
       #subtitle =,
       colour = "black") + 
  xlab("X coordinates (pixels)") + 
  ylab("Y coordinates (pixels)") + 
  theme(axis.title = element_text(size = rel(2)),
        axis.line.y.left = element_line(colour = "black"),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = rel(2)),
        axis.ticks.length.y.left = unit(.15, "cm"),
        axis.ticks.length.x.bottom = unit(.15, "cm"),
        plot.title = element_text(colour = "black", size = rel(2.5), hjust = 0.5),
        plot.subtitle = element_text(colour = "black", size = rel(1.9), hjust = 0.5),
        plot.margin = unit(c(2,3,2,3), "cm"),
        legend.title = element_text(colour = "black", size = rel(2)),
        legend.text = element_text(colour = "black", size = rel(1.7)),
        panel.background = element_rect(fill = "white"))
  


# ggplot(xy_coord, aes(x = Spot_X, y = Spot_Y, label = X_Y)) +
#   geom_point(size = 10) + 
#   xlim(15, 25) +
#   ylim(30, 50) +
#   geom_label()
# 

ggsave("voronoi_tessellation_cut2.pdf",
       width = grDevices::dev.size(units = "px")[1]/96,
       height = grDevices::dev.size(units = "px")[2]/96,
       units = "in",
       dpi = 400)

#---------------------TEST STUF...------------------------------#
#---------------------------------------------------------------#
