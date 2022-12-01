#' A series of code that will lead to the development of the GW-kNN method of 
#' clustering.
#'


# Make observations matrix g x s ----
## select the 500 most variable genes
obs <- counts[select,] %>% 
    .[,nb_names] %>%
    t() %>% 
    as.data.frame()

# Make spatial weights matrix s x s ----
## set parameters ----
dist.Mat <- gw.dist(dp.locat = st_coordinates(polygons$geom_cntd), p = 2)
bw = (range(dist.Mat)[2])+4
kernel = "bisquare"
## Calculate W distance matrix ----
w <- gw.weight(vdist = dist.Mat, 
               bw = bw,
               kernel = kernel,
               adaptive = FALSE)
## Plot a weighted distance example ----
ggplot() +
    geom_line(aes(x = dist.Mat[,1], y = w[,1]), 
              colour = "black",
              linewidth = 1.5) + 
    geom_point(aes(x = dist.Mat[,1], y = w[,1]), 
               colour = "skyblue", 
               alpha = 0.4) + 
    labs(title = kernel,
         subtitle = paste0("bw = ", bw),
         caption = paste0(scaled, " VST data")) +
    ylab("weight") + 
    xlab("distance") +
    my_theme

# Make the array s x g x s ----
# ## Set dimensions
# loc.n <- nrow(obs) # number of locations (y-axis)
# var.n <- ncol(obs) # number of genes (x-axis)
# focus.n <- loc.n   # number of iterations (z-axis)
# 
# ## -------------------------------------------- ##
# ### TO DO: 
# #### add dimnames to the array:
# #### dimnames=list(row.names,column.names,matrix.names)
# obs.W <- array(data = 0, c(loc.n, var.n, focus.n)) #3D array: [no of rows, no of cols, no of arrays]
# 
# # Weight expression data
# ## select focus point
# focus <- 1
# ## weight 
# w.counts <- get.gw.counts(obs = obs, wdmat = w, focus = focus)
# ## store
# obs.W[,,focus] <- w.counts
obs.W <- get.gwCount.array(obs, w)

# Calculate distances ----
dist.W <- dist(obs.W[,,1], method = "euclidean", p = 2) %>% 
    as.matrix()

# Get nearest neighbours ----