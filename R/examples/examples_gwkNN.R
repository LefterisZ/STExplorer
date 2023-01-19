#' A series of code that will lead to the development of the GW-kNN method of 
#' clustering.
#'

library(graphlayouts)
library(igraph)
library(ggraph)

# 1. Create a gwKNN graphics output directory ----
gwknnDir <- file.path(graphDir, "gwknn")

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
dist.W <- get.gwDist.array(obs.W, method = "euclidean", p = 2)

# Get nearest neighbours ----
knn.W <- get.gwKNN.list(dist.W, k = 7)

# Get the graphs ----
graph.W <- get.gwGraph.array(knn.W)

# Get edge frequencies ----
edge.Frq <- get.edge.freq(graph.W)

# Get edge stats ----
edges <- get.edge.stats(edge.Frq)
    ## filter out low-weight edges ----
    sum.stats <- summary(edges$weight) # get summary stats
    sum(edges$weight < sum.stats["1st Qu."]) # how many are below the 1st Qu. point?
    edges <- edges %>% 
        arrange(desc(weight)) %>%
        group_by(from) %>%
        slice(1:6)

    
# Prepare graph for clustering ----
graph <- graph.data.frame(edges, directed = FALSE)

E(graph)$weight <- edges$weight
V(graph)
I(graph)

# Cluster with Louvain ----
clusters_Louv <- cluster_louvain(graph = graph,
                                 weights = E(graph)$weight,
                                 resolution = 0.75)
membership(clusters_Louv)
V(graph)$louv <- as.factor(membership(clusters_Louv))
order(V(graph)$louv)
col.No = length(clusters_Louv)
colour.values <- get.colours(col.No)


#layout <- layout_as_tree(graph)
l = "stress."
prefix.louvain <- paste0("Louv")
res = "0.75"
filt = "filt_6NNs"
## Plot graph on an automated layout ----
ggraph::ggraph(g = graph,
               layout = "auto") + 
    geom_edge_link(aes(colour = E(graph)$weight)) + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = colour.values) +
    scale_edge_width(range = c(0.1, 2))
ggsave(file.path(gwknnDir, paste0("gwKNN.clust.", prefix.louvain, ".", res, ".", l, filt, ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Plot clusters on the tissue map ----
clusters_Louv_df <- clusters_Louv$names %>%
    as.data.frame() %>%
    mutate(cluster = as.factor(clusters_Louv$membership)) %>%
    rename("Barcode" = ".")

map <- clusters_Louv_df  %>%
    left_join(polygons[,c("Barcode", "geom_pol")])

ggplot() + 
    geom_sf(data = map$geom_pol,
            aes(fill = map$cluster)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "GWKNN with Louvain clusters",
         caption = paste0(prefix.louvain, " res = ", res,
                          "\n", filt),
         fill = "Clusters") + 
    my_theme
ggsave(file.path(gwknnDir, paste0("gwKNN.clust.", prefix.louvain, ".", res, ".map.", filt, ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Cluster with Leiden ----
clusters_Leid <- cluster_leiden(graph, 
                                objective_function = "modularity",
                                n_iterations = 2,
                                resolution_parameter = 0.75)
membership(clusters_Leid)
V(graph)$leid <- as.character(membership(clusters_Leid))
col.No = length(clusters_Leid)
colour.values <- get.colours(col.No)

#layout <- layout_as_tree(graph)
l = "stress."
prefix.leiden <- paste0("Leid")
## Plot graph on an automated layout ----
ggraph::ggraph(g = graph,
               layout = "auto") + 
    geom_edge_link(aes(colour = E(graph)$weight)) + 
    geom_node_point(aes(colour = leid)) + 
    scale_colour_manual(values = colour.values) +
    scale_edge_width(range = c(0.1, 2))
ggsave(file.path(gwknnDir, paste0("gwKNN.clust.", prefix.leiden, ".", l, ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Plot clusters on the tissue map ----
clusters_Leid_df <- clusters_Leid$names %>%
    as.data.frame() %>%
    mutate(cluster = as.factor(clusters_Leid$membership)) %>%
    rename("Barcode" = ".")

map <- clusters_Leid_df  %>%
    left_join(polygons[,c("Barcode", "geom_pol")])

ggplot() + 
    geom_sf(data = map$geom_pol,
            aes(fill = map$cluster)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "GWKNN with Leiden clusters",
         fill = "Clusters") + 
    my_theme
ggsave(file.path(gwknnDir, paste0("gwKNN.clust.", prefix.louvain, ".", res, ".map.", filt, ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
