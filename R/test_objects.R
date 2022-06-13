## test data.frame of random discrete coordinates ----
## out of random discrete numbers between 1 and 1000
test_df <- as.data.frame(tibble(X = c(rdunif(100, b=1000, a=1)), 
                                Y = c(rdunif(100, b=1000, a=1))))


## test tables of spot positions on a slide ----
## --> 10X Visium slide
# slide of: 5000 spots
x <- 0:99
y1 <- seq(0, 98, 2)
y2 <- seq(1, 99, 2)

test_tissue_5K <- tibble(x = rep(x, each=50), 
                      y = rep(c(y1, y2), 50))

# slide of: 10000 spots
x <- 0:199
y1 <- seq(0, 198, 2)
y2 <- seq(1, 199, 2)

test_tissue_10K <- tibble(x = rep(x, each=100), y = rep(c(y1, y2), 100))

## test tissue sample ----
test_tissue_sf <- test_tissue_5K %>% 
    st_as_sf(coords = c("x", "y"))

test_tissue_sample <- test_tissue_sf %>% 
    st_sample(size = 15,
              type = "random",
              exact = TRUE)

## test centroids file ----
test_centroids <- mob_centroids
head(test_centroids)


## test directory path ----
test_sDir <- "/home/b9047753/Documents/projects/Visium_MOB/spaceranger_outs/Olfactory_Bulb/Olfactory_Bulb_A1_Results/spatial"
test_scaleFactors <- "scalefactors_json.json"

## test 10X Visium tissue_position_list.csv fille ----
test_tissue_positions <- mob_input

