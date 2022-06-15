# This script will deal with extracting information needed for neighbour
#   identification and weighted matrix calculation from the Spaniel sce
#   object.

# packages needed
library(spdep)
library(sf)
library(jsonlite)
library(ggplot2)

# non-package functions needed
source("./R/spot_diameter.R")
source("./R/make_bb_polygon.R")
source("./R/readSpaceranger.R")
source("./R/add_perimeter.R")
source("./R/spot_neighbours.R")

## set the file path to spaceranger's spatial folder
spatialDir = file.path(inputDir, "Olfactory_Bulb/Olfactory_Bulb_A1_Results/spatial")

## Import the dataset
input <- readSpaceranger(spatialDir, res = "low") %>% #read-in data
    add_perimeter() #add perimeter of spots around the tissue 

## Select spots in bins 1 and 2
spot_position <- input %>% 
    filter(new_bin == 1 | new_bin == 2) %>% 
    select(c("Barcode", "pixel_x", "pixel_y", "new_bin")) %>% 
    remove_rownames() 

## Convert spots to polygon centroids
centroids <- spot_position %>% 
  st_as_sf(coords = c("pixel_x", "pixel_y"))


## Combine the points into a multipoint geometry:
cntd_union <- st_union(centroids)
head(cntd_union)

## Use the union of points generate a voronoi object
voronoi <- st_voronoi(cntd_union, bOnlyEdges = TRUE)
head(voronoi)

## Create an enveloped voronoi tessellation around the tissue
voronoi_env <- st_intersection(st_cast(voronoi), st_convex_hull(cntd_union))

## plot the voronoi tessellation
ggplot() +
  geom_sf(data = voronoi_env, colour = "black", fill = "white") + 
  geom_sf(data = centroids, colour = centroids$new_bin) + 
  #xlim(boxXmin, boxXmax) + 
  #ylim(boxYmin, boxYmax) + 
  # Add titles and visually format the plot:
  labs(title = paste("Voronoi tessellation"),
       subtitle = ,
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


ggsave("voronoi_tessellation_with_perimeter.pdf",
       width = grDevices::dev.size(units = "px")[1]/96,
       height = grDevices::dev.size(units = "px")[2]/96,
       units = "in",
       dpi = 400)

#---------------------TEST STUF...------------------------------#
#---------------------------------------------------------------#
