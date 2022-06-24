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
source("./R/get_nb_IDs.R")

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
  st_as_sf(coords = c("pixel_x", "pixel_y"), 
           remove = FALSE)

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

## Generate the POLYGONS from the MULTILINESTRING and attach the barcode names 
##  from bin_1 only
polygons <- st_polygonize(voronoi_env) %>% # polygonise the tessellation
    st_cast() %>% # convert GEOMETRYCOLLECTION to multiple POLYGONS
    st_sf() %>%  # convert sfc object to sf for st_join afterwards
    st_join(., 
            centroids[centroids$new_bin == 1,],
            join = st_contains,
            left = FALSE) %>% # Join the centroids with the POLYGONS
    mutate(Barcode_rn = Barcode) %>% # duplicate the barcode column
    column_to_rownames("Barcode_rn") %>% # move duplicate column to row names
    st_sf() # convert back to sf (mutate makes it a df)

## Create contiguity neighbours
neighbours <- poly2nb(polygons, snap = 0)
names(neighbours) = attr(neighbours, "region.id") # add names to the sublists

## Add number of neighbours for each polygon back to the polygons object
polygons$nb_count <- card(neighbours)

## Add the neighbour IDs as a nested df in the polygons object
nb_IDs <- nb_sf %>%
    st_drop_geometry() %>%
    select(i_ID, j_ID) %>%
    rename(nb_IDs = j_ID) %>%
    group_by(i_ID) %>%
    nest()

polygons <- right_join(polygons, nbs, by = c("Barcode" = "i_ID")) %>%
    rename(nb_IDs = data)

## Get a neighbours object for ggplot2 plotting
spot_coord <- polygons[,c("pixel_x", "pixel_y")] %>%
    st_drop_geometry() # get the spot coords as a df

nb_sf <- as(nb2lines(neighbours, coords = spot_coord), "sf")

## Plot neighbours graph
ggplot() +
    geom_sf(data = polygons$geometry, colour = "grey30", fill = "white") +
    geom_sf(data = nb_sf, colour = "black") + 
    geom_point(data = polygons, aes(x = pixel_x, y = pixel_y, colour = factor(nb_count))) + 
    # Add titles and visually format the plot:
    scale_fill_discrete(c("#34568B", "#FF6F61", "#6B5B95", 
                          "#88B04B", "#F7CAC9", "#92A8D1")) +
    labs(title = paste("Contiguity neighbours"),
         subtitle = "",
         colour = "Neighbours\ncount") + 
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

ggsave("neighbours_graph-contiguity.pdf",
       width = grDevices::dev.size(units = "px")[1]/96,
       height = grDevices::dev.size(units = "px")[2]/96,
       units = "in",
       dpi = 400)

#---------------------TEST STUF...------------------------------#
#---------------------------------------------------------------#

