# In this script we try to simulate random tissue shapes to test them




#---------------------------------------------------------
# Generate a grid of spots resembling the 10X Visium slide
#---------------------------------------------------------
simulatissue <- make_slide(spots_x = 100, spots_y = 100)


## plot the spots
ggplot(tissue, aes(x=x, y=y)) + geom_point()

## convert the spots to an sf object
tissue_sf <- tissue %>% 
  st_as_sf(coords = c("x", "y"))

## take a random set from the spots and sort them based on a theoretical circle
tissue_sample <- tissue_sf %>% 
    st_sample(size = 15, 
              type = "random", 
              exact = TRUE) %>% 
    spot_sort() %>% 
    spot_add(return = "sfc")

##-----------------------------------------------##
## a good set of spots for later                 ##
## ddf = df                                      ##
##-----------------------------------------------##

ggplot() + 
    geom_sf(data = tissue_sf, colour = "black") + 
    geom_sf(data = tissue_sample, colour = "red2", size = 2.5) + 
    geom_sf(data = st_cast(tissue_sample, "MULTILINESTRING"), 
            colour = "dark blue",
            size = 1)

ggsave(file.path(graphDir, "random_tissue_shape.pdf"),
       width = grDevices::dev.size(units = "px")[1]/96,
       height = grDevices::dev.size(units = "px")[2]/96,
       units = "in",
       dpi = 600)

