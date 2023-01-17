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
obs.W <- get.gwCount.array(obs, w)

# Calculate distances ----
dist.W <- get.gwDist.array(obs.W, focus.n = c(1, 2, 3), method = "euclidean", p = 2)

# Get nearest neighbours ----
knn.W <- get.gwKNN.list(dist.W, k = 7, focus.n = c(1, 2, 3))

# Get the graphs ----
graph.W <- get.gwGraph.array(knn.W, focus.n = c(1, 2, 3))

# Get edge frequencies ----
edge.Frq <- get.edge.freq(graph.W, k = 7)


