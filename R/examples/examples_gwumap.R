install.packages("umap")
install.packages("uwot")
install.packages("dbscan")
install.packages("graphlayouts")
install.packages("ggraph")
install.packages("snahelper")
library(umap)
library(uwot)
library(dbscan)
library(graphlayouts)
library(igraph)
library(ggraph)


# 1. Create a gwumap graphics output directory ----
gwumapDir <- file.path(graphDir, "gwumap")

# 2a. Create input data from 500 most variable genes ----
inputUMAP <- vst_df[, select] %>%
    as.data.frame() %>% 
    .[nb_names,]
# 2b. Create input data from the first 50 Principal Components
## princomp function was used for PCA on 500 most variable genes
### input to princomp was scaled
inputUMAP.pca <- pca[["scores"]][,1:30] %>% 
    as.data.frame() %>%
    .[nb_names,]

# ---------------------------------------------------------------------------- #
# 3. Run UMAP - default - Un-scaled VST data ----
umap <- umap(inputUMAP) # generate UMAP

umap.layout_df <- as.data.frame(umap$layout)
head(umap.layout_df)
polygons[1:6, 1:2]

ggplot() + 
    geom_point(data = umap.layout_df,
               aes(x = V1, y = V2)) +
    labs(title = "Un-scaled VST data UMAP") + 
    my_theme

prefix.gwumap <- "gwumap_default_500.vst.unscaled"
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# ---------------------------------------------------------------------------- #
# 4. Run UMAP - default - Scaled VST data ----
inputUMAP.scaled <- scale(inputUMAP) # scale data by column (gene)
umap.scaled <- umap(inputUMAP.scaled) # generate UMAP

umap.layout_df.scaled <- as.data.frame(umap.scaled$layout)
head(umap.layout_df.scaled)
polygons[1:6, 1:2]

ggplot() + 
    geom_point(data = umap.layout_df.scaled,
               aes(x = V1, y = V2)) +
    labs(title = "Scaled VST data UMAP") + 
    my_theme

prefix.gwumap <- "gwumap_default_500.vst.scaled"
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# ---------------------------------------------------------------------------- #
# 5. Run UMAP - custom.1 - Scaled VST data ----
## Customise configuration ----
custom.config.1 <- umap.defaults
custom.config.1$n_neighbors <- 6
## Run UMAP ----
umap.custom.1 <- umap(inputUMAP.scaled,
                      config = custom.config.1)

umap.layout_df.1 <- as.data.frame(umap.custom.1$layout)
head(umap.layout_df.1)
polygons[1:6, 1:2]
## Plot UMAP ----
ggplot() + 
    geom_point(data = umap.layout_df.1,
               aes(x = V1, y = V2)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 6") + 
    my_theme

prefix.gwumap <- "gwumap_cust.n6_500.vst.unscaled"
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
## Manually group, plot and map ----
umap.map.1 <- umap.layout_df.1 %>%
    mutate("group" = ifelse(.data$V2 < 2.2, "1", 
                            ifelse(.data$V2 > 5.5 & .data$V1 < 0, "2", 
                                   ifelse(.data$V2 > 5.5 & .data$V1 >0, "3", "4")))) %>%
    mutate("geometry" = polygons$geom_pol)

ggplot() + 
    geom_point(data = umap.map.1,
               aes(x = V1, y = V2, colour = group)) + 
    scale_color_manual(values = c4a("wright25", 4)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 6") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.grouped.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

ggplot() + 
    geom_sf(data = umap.map.1$geometry,
            aes(fill = umap.map.1$group)) + 
    scale_fill_manual(values = c4a("wright25", 4)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 6") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.map.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# ---------------------------------------------------------------------------- #
# 6. Run UMAP - custom.2 - Scaled VST data ----
## Customise configuration ----
custom.config.2 <- umap.defaults
custom.config.2$n_neighbors <- 6
custom.config.2$min_dist <- 0.05
## Run UMAP ----
umap.custom.2 <- umap(inputUMAP.scaled,
                      config = custom.config.2)

umap.layout_df.2 <- as.data.frame(umap.custom.2$layout)
head(umap.layout_df.2)
polygons[1:6, 1:2]
## Plot UMAP ----
ggplot() + 
    geom_point(data = umap.layout_df.2,
               aes(x = V1, y = V2)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 6\nmin dist = 0.05") + 
    my_theme

prefix.gwumap <- "gwumap_cust.n6.dist005_500.vst.unscaled"
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
## Manually group, plot and map ----
umap.map.2 <- umap.layout_df.2 %>%
    mutate("group" = ifelse(.data$V2 < 0, "1", 
                            ifelse(.data$V2 > 8, "2", 
                                   ifelse(.data$V2 > 0 & .data$V2 < 7.5 & .data$V1 > 1, "3", "4")))) %>%
    mutate("geometry" = polygons$geom_pol)

ggplot() + 
    geom_point(data = umap.map.2,
               aes(x = V1, y = V2, colour = group)) + 
    scale_color_manual(values = c4a("wright25", 4)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 6\nmin dist = 0.05") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.grouped.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

ggplot() + 
    geom_sf(data = umap.map.2$geometry,
            aes(fill = umap.map.2$group)) + 
    scale_fill_manual(values = c4a("wright25", 4)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 6\nmin dist = 0.05") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.map.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# ---------------------------------------------------------------------------- #
# 7. Run UMAP - custom.3 - Scaled VST data ----
## Customise configuration ----
custom.config.3 <- umap.defaults
custom.config.3$n_neighbors <- 7
custom.config.3$min_dist <- 0.05
## Run UMAP ----
umap.custom.3 <- umap(inputUMAP.scaled,
                      config = custom.config.3)

umap.layout_df.3 <- as.data.frame(umap.custom.3$layout)
head(umap.layout_df.3)
polygons[1:6, 1:2]
## Plot UMAP ----
ggplot() + 
    geom_point(data = umap.layout_df.3,
               aes(x = V1, y = V2)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 7\nmin dist = 0.05") + 
    my_theme

prefix.gwumap <- "gwumap_cust.n7.dist005_500.vst.unscaled"
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
## Manually group, plot and map ----
umap.map.3 <- umap.layout_df.3 %>%
    mutate("group" = ifelse(.data$V2 > 0 & .data$V2 < 2.5 & .data$V1 < 1.8 & .data$V1 > 0.5, "5",
                            ifelse(.data$V2 > 0, "1", 
                                   ifelse(.data$V2 < 0 & .data$V2 > -6, "3", 
                                          ifelse(.data$V2 < -8.2 & .data$V1 > -0.8 & .data$V1 < 0.4, "2", "4"))))) %>%
    mutate("geometry" = polygons$geom_pol)

ggplot() + 
    geom_point(data = umap.map.3,
               aes(x = V1, y = V2, colour = group)) + 
    scale_color_manual(values = c4a("wright25", 7)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 7\nmin dist = 0.05") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.grouped.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

ggplot() + 
    geom_sf(data = umap.map.3$geometry,
            aes(fill = umap.map.3$group)) + 
    scale_fill_manual(values = c4a("wright25", 7)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 7\nmin dist = 0.05") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.map.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## kmeans group, plot and map ----
kmeans.7 <- kmeans(umap.layout_df.3, centers = 7)
head(umap.layout_df.3)
kmeans.7$cluster[1:6]
umap.map.3.kmean <- umap.layout_df.3 %>%
    mutate("Cluster" = as.factor(kmeans.7$cluster)) %>%
    mutate("geometry" = polygons$geom_pol)

ggplot() + 
    geom_point(data = umap.map.3.kmean,
               aes(x = V1, y = V2, colour = Cluster)) + 
    scale_color_manual(values = c4a("wright25", 7)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 7\nmin dist = 0.05\nkmeans = 7") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.grouped.kmean.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

ggplot() + 
    geom_sf(data = umap.map.3.kmean$geometry,
            aes(fill = umap.map.3.kmean$Cluster)) + 
    scale_fill_manual(values = c4a("wright25", 7)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 7\nmin dist = 0.05\nkmeans = 7") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.map.kmean.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Plot with igraph ----
### Prepare a from-to data.frame ----
umap.map.3.igraph <- umap.custom.3$knn$indexes
colnames(umap.map.3.igraph) <- c("from", 
                                 paste0("nb", 1:(dim(umap.map.3.igraph)[2]-1)))
umap.map.3.igraph <- umap.map.3.igraph %>%
    as.data.frame() %>%
    pivot_longer(-from, names_to = NULL, values_to = "to")

### Make the igraph object ----
umap.map.3.ig <- graph.data.frame(umap.map.3.igraph, directed = FALSE)

E(umap.map.3.ig)
V(umap.map.3.ig)$grp <- umap.map.3$group
I(umap.map.3.ig)

ggraph::ggraph(g = umap.map.3.ig,
               layout = "manual",
               x = umap.custom.3$layout[,1],
               y = umap.custom.3$layout[,2]) + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = grp)) + 
    scale_colour_manual(values = c4a("wright25", 7))
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.igraph.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

ggraph::ggraph(g = umap.map.3.ig,
               layout = "stress") + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = grp)) + 
    scale_colour_manual(values = c4a("wright25", 7))
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.igraph2.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

### Cluster with Louvain ----
umap.map.3.clst <- cluster_louvain(umap.map.3.ig)
membership(umap.map.3.clst)
V(umap.map.3.ig)$louv <- as.character(membership(umap.map.3.clst))

ggraph::ggraph(g = umap.map.3.ig,
               layout = "manual",
               x = umap.custom.3$layout[,1],
               y = umap.custom.3$layout[,2]) + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = c4a("wright25", 9))
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.igraph.louv.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

ggraph::ggraph(g = umap.map.3.ig,
               layout = "stress") + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = c4a("wright25", 9))
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.igraph.louv2.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

umap.map.3 <- umap.map.3 %>%
    mutate("louvain" = as.factor(membership(umap.map.3.clst)))

ggplot() + 
    geom_sf(data = umap.map.3$geometry,
            aes(fill = umap.map.3$louvain)) + 
    scale_fill_manual(values = c4a("wright25", 9)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 7\nmin dist = 0.05",
         fill = "Louvain clusters") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.map.louvain.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

### Apply spatial weights to graph edges
dist.Mat <- gw.dist(dp.locat = st_coordinates(polygons$geom_cntd),
                    p = 2)
a = 1
b = 1
dist.Mat.w1 <- ((1/(dist.Mat+a))+b)

bw = 3*spot_diameter(spatialDir)
kernel = "bisquare"

dist.Mat.w2 <- gw.weight(vdist = dist.Mat, 
                            bw = bw,
                            kernel = kernel,
                            adaptive = FALSE)

umap.map.3.edgeW1 <- apply(umap.map.3.igraph, 1, function(x){
    dist.Mat.w1[x[1],x[2]]
})
umap.map.3.edgeW2 <- apply(umap.map.3.igraph, 1, function(x){
    dist.Mat.w2[x[1],x[2]]
})

range(umap.map.3.edgeW1)
E(umap.map.3.ig)$weights <- umap.map.3.edgeW1


umap.map.3.ig.filt <- subgraph.edges(umap.map.3.ig, 
                                     E(umap.map.3.ig)[E(umap.map.3.ig)$weights > 1], 
                                     del=F)

umap.map.3.clst.filt <- cluster_louvain(umap.map.3.ig.filt)
membership(umap.map.3.clst.filt)
V(umap.map.3.ig.filt)$louv <- as.factor(membership(umap.map.3.clst.filt))

ggraph::ggraph(g = umap.map.3.ig.filt,
               layout = "manual",
               x = umap.custom.3$layout[,1],
               y = umap.custom.3$layout[,2]) + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = c4a("wright25", 9))

## Apply spatial weights to UMAP layout ----
### Make sf point object ----
umap.map.3.cntd <- umap.custom.3$layout %>% # retrieve the UMAP layout
    as.data.frame() %>%
    mutate(Barcode = rownames(umap.custom.3$layout))
colnames(umap.map.3.cntd) <- c("UMAP1", "UMAP2", "Barcode") # add colnames

umap.map.3.cntd <- umap.map.3.cntd %>% # make the sf point object
    st_as_sf(coords = c("UMAP1", "UMAP2"), 
             remove = FALSE)

umap.map.3.cntd[1:6,]
polygons[1:6,1:3]

### Find distances between points  ----
umap.map.3.dist <- gw.dist(dp.locat = st_coordinates(umap.map.3.cntd),
                           p = 2)
range(umap.map.3.dist)

### Calculate weighted distances ----
dist.Mat <- gw.dist(dp.locat = st_coordinates(polygons$geom_cntd),
                    p = 2)
bw = 3*spot_diameter(spatialDir)
kernel = "exponential"

polygons.wdist <- gw.weight(vdist = dist.Mat, 
                            bw = bw,
                            kernel = kernel,
                            adaptive = FALSE)
range(polygons.wdist)
range(dist.Mat)

### Calculate weighted UMAP distances
umap.map.3.wdist <- umap.map.3.dist / polygons.wdist

umap.map.3.cmdscale <- cmdscale(umap.map.3.wdist, k = 2, add = TRUE) %>%
    .[["points"]] %>%
    as.data.frame() %>%
    mutate(group = umap.map.3$group) %>%
    rename("UMAP1" = "V1", "UMAP2" = "V2")

### Plot the weighted UMAP layout ----
ggplot() + 
    geom_point(data = umap.map.3.cmdscale,
               aes(x = UMAP1, y = UMAP2, colour = group)) + 
    scale_color_manual(values = c4a("wright25", 7)) +
    labs(title = "Geographically weighted UMAP",
         subtitle = "n neighbours = 7\nmin dist = 0.05") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.gwumap.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

### Reconstruct UMAP from point distances ----
cmdscale <- cmdscale(umap.map.3.dist, k = 2, add = TRUE) %>%
    .[["points"]] %>%
    as.data.frame() %>%
    mutate(group = umap.map.3$group) %>%
    rename("UMAP1" = "V1", "UMAP2" = "V2")
ggplot() +
    geom_point(data = cmdscale,
               aes(x = UMAP1, y = UMAP2, colour = group)) +
    scale_color_manual(values = c4a("wright25", 7)) +
    labs(title = "Reconstructed UMAP from point distances",
         subtitle = "n neighbours = 7\nmin dist = 0.05") +
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.reconst.umap.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
### Reconstruct tissue map from weighted point distances ----
cmdscale_polygons.w <- cmdscale(polygons.wdist, k = 2, add = TRUE) %>%
    .[["points"]] %>%
    as.data.frame() %>%
    mutate(group = umap.map.3$group) %>%
    rename("X" = "V1", "Y" = "V2")
ggplot() +
    geom_point(data = cmdscale_polygons.w,
               aes(x = X, y = Y, colour = group)) +
    scale_color_manual(values = c4a("wright25", 7)) +
    labs(title = "Reconstructed map from \nweighted point distances",
         subtitle = "kernel = exponential") +
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.reconst.mapW.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

### Reconstruct tissue map from simple point distances ----
cmdscale_polygons <- cmdscale(dist.Mat, k = 2, add = TRUE) %>%
    .[["points"]] %>%
    as.data.frame() %>%
    mutate(group = umap.map.3$group) %>%
    rename("X" = "V1", "Y" = "V2")
ggplot() +
    geom_point(data = cmdscale_polygons,
               aes(x = X, y = Y, colour = group), size = 5.8) +
    scale_color_manual(values = c4a("wright25", 7)) +
    labs(title = "Reconstructed map from point distances") +
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.reconst.map.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Calculate knn from UMAP layout ----
umap.custom.3.knn <- knearneigh(umap.custom.3$layout, k = 6)

# ---------------------------------------------------------------------------- #
# 8. Run UMAP - custom.4 - Scaled VST data ----
## Customise configuration ----
custom.config.4 <- umap.defaults
custom.config.4$n_neighbors <- 50
custom.config.4$min_dist <- 0.05
## Run UMAP ----
umap.custom.4 <- umap(inputUMAP.scaled,
                      config = custom.config.4)

umap.layout_df.4 <- as.data.frame(umap.custom.4$layout)
head(umap.layout_df.4)
polygons[1:6, 1:2]
## Plot UMAP ----
ggplot() + 
    geom_point(data = umap.layout_df.4,
               aes(x = V1, y = V2)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = "n neighbours = 50\nmin dist = 0.05") + 
    my_theme

prefix.gwumap <- "gwumap_cust.n7.dist005_500.vst.unscaled"
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# ---------------------------------------------------------------------------- #
# 9. Run UMAP - neighbours = all spots ----
## Customise configuration ----
custom.config.5 <- umap.defaults
custom.config.5$n_neighbors <- dim(inputUMAP.scaled)[1]
custom.config.5$min_dist <- 0.01

## Run UMAP ----
scaled <- "scaled"
umap.custom.5 <- umap(inputUMAP.scaled,
                    config = custom.config.5)

umap.layout_df.5 <- as.data.frame(umap.custom.5$layout)
head(umap.layout_df.5)
polygons[1:6, 1:2]

## Plot UMAP ----
ggplot() + 
    geom_point(data = umap.layout_df.5,
               aes(x = V1, y = V2)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = paste0("knn = ", custom.config.5$n_neighbors,
                           "\nminDist = ", custom.config.5$min_dist)) + 
    my_theme

prefix.gwumap <- paste0("gwumap_cust.n", custom.config.5$n_neighbors,
                        ".dist", custom.config.5$min_dist, ".", scaled)
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Get the knn indexes and the UMAP distances ----
umap.map.5.igraph.idx <- umap.custom.5$knn$indexes
colnames(umap.map.5.igraph.idx) <- c("from", 
                                   paste0("nb", 1:(dim(umap.map.5.igraph.idx)[2]-1)))

umap.map.5.igraph.dist <- umap.custom.5$knn$distances
colnames(umap.map.5.igraph.dist) <- c("from", 
                                    paste0("nb", 1:(dim(umap.map.5.igraph.dist)[2]-1)))

## Build EDGES from-to data.frame ----
umap.map.5.igraph.E <- umap.map.5.igraph.idx %>%
    as.data.frame() %>%
    pivot_longer(-from, names_to = NULL, values_to = "to")

## Plot UN-weighted LOUVAIN ----
### Generate the graph object ----
umap.map.5.ig <- graph.data.frame(umap.map.5.igraph.E, directed = FALSE)
umap.map.5.clst <- cluster_louvain(umap.map.5.ig)
V(umap.map.5.ig)$louv <- as.character(membership(umap.map.5.clst))
### Plot ----
ggraph::ggraph(g = umap.map.5.ig,
               layout = "manual",
               x = umap.custom.5$layout[,1],
               y = umap.custom.5$layout[,2]) + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = c4a("wright25", 11))
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, ".louv.default.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

umap.map.5.igraph.D <- umap.map.5.igraph.dist %>%
    as.data.frame() %>%
    pivot_longer(-from, names_to = NULL, values_to = "dist") %>%
    as.data.frame() %>%
    select(-from) %>%
    .[,1]
E(umap.map.5.ig)$weight <- umap.map.5.igraph.D
umap.map.5.clst <- cluster_louvain(umap.map.5.ig)
V(umap.map.5.ig)$louv <- as.character(membership(umap.map.5.clst))


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# 10. Run GW-UMAP - weight knn statistical distance ----
## Customise configuration ----
custom.config <- umap.defaults
custom.config$n_neighbors <- 10
custom.config$min_dist <- 0.2
custom.config$n_components <- 2

## Run UMAP ----
scaled <- "scaled"
umap.custom <- umap::umap(inputUMAP.scaled,
                      config = custom.config)

umap.layout_df <- as.data.frame(umap.custom$layout)
head(umap.layout_df)
polygons[1:6, 1:2]

## Plot UMAP ----
ggplot() + 
    geom_point(data = umap.layout_df,
               aes(x = V1, y = V2)) +
    labs(title = "Scaled VST data UMAP",
         subtitle = paste0("knn = ", custom.config$n_neighbors,
                           "\nminDist = ", custom.config$min_dist)) + 
    my_theme

prefix.gwumap <- paste0("gwumap_cust.n", custom.config$n_neighbors,
                        ".dist", custom.config$min_dist, ".", scaled)
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, "_embedding.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Get the knn indexes and the UMAP distances ----
umap.map.igraph.idx <- umap.custom$knn$indexes
colnames(umap.map.igraph.idx) <- c("from", 
                                     paste0("nb", 1:(dim(umap.map.igraph.idx)[2]-1)))

umap.map.igraph.dist <- umap.custom$knn$distances
colnames(umap.map.igraph.dist) <- c("from", 
                                      paste0("nb", 1:(dim(umap.map.igraph.dist)[2]-1)))

## Build EDGES from-to data.frame ----
umap.map.igraph.E <- umap.map.igraph.idx %>%
    as.data.frame() %>%
    pivot_longer(-from, names_to = NULL, values_to = "to")

## Plot UN-weighted LOUVAIN ----
### Generate the graph object ----
umap.map.ig <- graph.data.frame(umap.map.igraph.E, directed = FALSE)
umap.map.clst <- cluster_louvain(umap.map.ig)
V(umap.map.ig)$louv <- as.character(membership(umap.map.clst))
### Plot UMAP layout----
ggraph::ggraph(g = umap.map.ig,
               layout = "manual",
               x = umap.custom$layout[,1],
               y = umap.custom$layout[,2]) + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = c4a("wright25", 11))
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, ".louv.default.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
### Plot map layout----
umap.map <- umap.layout_df  %>%
    mutate("geometry" = polygons$geom_pol) %>%
    mutate("louvain" = as.factor(membership(umap.map.clst)))

ggplot() + 
    geom_sf(data = umap.map$geometry,
            aes(fill = umap.map$louvain)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "Louvain clusters",
         caption = paste0("knn = ", custom.config$n_neighbors,
                          "\nminDist = ", custom.config$min_dist,
                          "\nscaled = ", scaled),
         fill = "Louvain clusters") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.gwumap, ".map.louv.default.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Build EDGES UMAP distance vector ----
umap.map.igraph.D <- umap.map.igraph.dist %>%
    as.data.frame() %>%
    pivot_longer(-from, names_to = NULL, values_to = "dist") %>%
    as.data.frame() %>%
    select(-from) %>%
    .[,1]

## Calculate weighted physical distance matrix ----
### Set parameters ----
dist.Mat <- gw.dist(dp.locat = st_coordinates(polygons$geom_cntd), p = 2)
# a = 1
# b = 1
bw = (range(dist.Mat)[2])+4
kernel = "bisquare"
### Calculate W distance matrix ----
# dist.Mat.w1 <- 1/dist.Mat
dist.Mat.w2 <- gw.weight(vdist = dist.Mat, 
                         bw = bw,
                         kernel = kernel,
                         adaptive = FALSE)
### Plot weights distribution for a single location ----
ggplot() +
    geom_line(aes(x = dist.Mat[,1], y = dist.Mat.w2[,1]), 
              colour = "black",
              size = 1.5) + 
    geom_point(aes(x = dist.Mat[,1], y = dist.Mat.w2[,1]), 
               colour = "skyblue", 
               alpha = 0.4) + 
    labs(title = kernel,
         subtitle = paste0("bw = ", bw),
         caption = paste0(scaled, " VST data")) +
    ylab("weight") + 
    xlab("distance") +
    my_theme
ggsave(file.path(gwumapDir, 
                 paste0("dist.weight.curve.", scaled, ".", kernel, ".bw", round(bw), ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Build WEIGHTS vector for the graph edges ----
umap.map.igraph.W <- get.edge.weights(umap.map.igraph.E, dist.Mat.w2)
range(umap.map.igraph.W)

## Calculate WEIGHTED UMAP DISTANCES vector ----
umap.map.igraph.DW <- umap.map.igraph.D * umap.map.igraph.W
range(umap.map.igraph.DW)

## Get some stats about distances ----
umap.map.igraph.stats <- data.frame(umap.map.igraph.D, 
                                    umap.map.igraph.W,
                                    umap.map.igraph.DW)
summary(umap.map.igraph.stats)

## Add weighted UMAP distances as edge weights ----
umap.map.ig <- graph.data.frame(umap.map.igraph.E, directed = FALSE)
E(umap.map.ig)$weight <- umap.map.igraph.DW
V(umap.map.ig)
I(umap.map.ig)

## LOUVAIN CLUSTERING with EDGE WEIGHTS ----
umap.map.clst.DW <- cluster_louvain(umap.map.ig, weights = E(umap.map.ig)$weight)
membership(umap.map.clst.DW)
V(umap.map.ig)$louv <- as.character(membership(umap.map.clst.DW))
order(V(umap.map.ig)$louv)
col.No = length(umap.map.clst.DW)
colour.values <- get.colours(col.No)

### calculate a layout with weights ----
layout <- layout_as_tree(umap.map.ig)
l = ""
### Set graph output prefix ----
prefix.louvain <- paste0("umap.custom.louv_embedding.DW.")

### Plot graph on the UMAP layout ----
ggraph::ggraph(g = umap.map.ig,
               layout = "auto",
                # x = umap.custom$layout[,1],
                # y = umap.custom$layout[,2]
               ) + 
    geom_edge_link(aes(#width = E(umap.map.ig)$weight,
                       colour = E(umap.map.ig)$weight), 
                   # edge_colour = "grey66"
                   ) + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = colour.values) +
    scale_edge_width(range = c(0.1, 2))
ggsave(file.path(gwumapDir, paste0(prefix.louvain, scaled, ".", kernel, ".bw", round(bw), ".", l, ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

### Plot graph on the tissue map ----
umap.map <- umap.layout_df  %>%
    mutate("geometry" = polygons$geom_pol) %>%
    mutate("louvain" = as.factor(membership(umap.map.clst.DW)))

ggplot() + 
    geom_sf(data = umap.map$geometry,
            aes(fill = umap.map$louvain)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "Louvain clusters",
         caption = paste0("knn = ", custom.config$n_neighbors,
                           "\nminDist = ", custom.config$min_dist,
                           "\nkernel = ", kernel,
                           "\nbw = ", round(bw),
                           "\nscaled = ", scaled),
         fill = "Louvain clusters") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.louvain, ".map.", scaled, ".", kernel, ".bw", round(bw), ".", l, ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## LEIDEN CLUSTERING with EDGE WEIGHTS ----
umap.map.clst.DW <- cluster_leiden(umap.map.ig, 
                                   objective_function = "modularity",
                                   n_iterations = 2,
                                   resolution_parameter = 1)
membership(umap.map.clst.DW)
V(umap.map.ig)$leid <- as.character(membership(umap.map.clst.DW))
col.No = length(umap.map.clst.DW)
colour.values <- get.colours(col.No)
### Set graph output prefix ----
prefix.leiden <- paste0("umap.custom.leid_embedding.DW.")

### Plot graph on the UMAP layout ----
ggraph::ggraph(g = umap.map.ig,
               layout = "manual",
               x = umap.custom$layout[,1],
               y = umap.custom$layout[,2]) + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = leid)) + 
    scale_colour_manual(values = colour.values)
ggsave(file.path(gwumapDir, paste0(prefix.leiden, scaled, ".", kernel, ".bw", round(bw), ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

### Plot graph on the tissue map ----
umap.map <- umap.layout_df %>%
    mutate("geometry" = polygons$geom_pol) %>%
    mutate("leiden" = as.factor(membership(umap.map.clst.DW)))

ggplot() + 
    geom_sf(data = umap.map$geometry,
            aes(fill = umap.map$leiden)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "Leiden clusters",
         caption = paste0("knn = ", custom.config$n_neighbors,
                          "\nminDist = ", custom.config$min_dist,
                          "\nkernel = ", kernel,
                          "\nbw = ", round(bw),
                          "\nscaled = ", scaled),
         fill = "Leiden clusters") + 
    my_theme
ggsave(file.path(gwumapDir, paste0(prefix.leiden, ".map.", scaled, ".", kernel, ".bw", round(bw), ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## HDBSCAN CLUSTERING of UMAP embedding ----
## Customise configuration ----
data.in <- expand.grid(knn = c(7, 10, 13, 15, 20),
                       minDist = c(0.001, 0.01, 0.1, 0.25, 0.5),
                       lowDims = c(2, 3, 5, 10, 15),
                       minPTS = c(8, 10, 15),
                       stringsAsFactors = FALSE)

apply(data.in, 1, function(x){
    custom.config <- umap.defaults
    custom.config$n_neighbors <- x[1]
    custom.config$min_dist <- x[2]
    custom.config$n_components <- x[3]

    ## Run UMAP ----
    scaled <- "scaled"
    umap.custom <- umap::umap(inputUMAP.pca,
                              config = custom.config)
    
    ## Run HFBSCAN ----
    umap.hdbsc.in <- data.frame(umap.custom$layout[,])
    colnames(umap.hdbsc.in) <- paste0("UMAP", 1:dim(umap.hdbsc.in)[2])
    
    umap.map.clst.HDBSC <- hdbscan(umap.hdbsc.in, 
                                   minPts = x[4],
                                   verbose = TRUE)
    
    #minPts = umap.map.clst.HDBSC$hc$call$minPts
    minPts = x[4]
    knn = custom.config$n_neighbors
    minDist = custom.config$min_dist
    lowD = custom.config$n_components
    col.No = length(unique(umap.map.clst.HDBSC$cluster))
    colour.values <- get.colours(col.No)
    ### Set graph output prefix ----
    prefix.hdbscan <- paste0("umap.custom.HDBSCAN.")
    
    ### Plot graph on the tissue map ----
    umap.map <- umap.hdbsc.in %>%
        mutate("geometry" = polygons$geom_pol) %>%
        mutate( "hdbscan" = as.factor(umap.map.clst.HDBSC$cluster))
    
    p <- ggplot() + 
        geom_sf(data = umap.map$geometry,
                aes(fill = umap.map$hdbscan)) + 
        scale_fill_manual(values = colour.values) +
        labs(title = "HDBSCAN clusters",
             caption = paste0("knn = ", custom.config$n_neighbors,
                              "\nminDist = ", custom.config$min_dist,
                              "\nminPts = ", minPts,
                              "\nlowDim = ", lowD,
                              "\nscaled = ", scaled),
             fill = "HDBSCAN clusters") + 
        my_theme
    print(p)
    ggsave(file.path(gwumapDir, paste0(prefix.hdbscan, minPts, ".map.pca.", scaled, ".knn", knn, ".minD", minDist, ".lowDim", lowD, ".tiff")),
           plot = p,
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 150)
})

## FUZZY CLUSTERING of UMAP embedding ----
## Customise configuration ----
data.in <- expand.grid(knn = c(7, 10, 13, 15, 20),
                       minDist = c(0.001, 0.01, 0.1, 0.25, 0.5),
                       lowDims = c(2, 3, 5, 10, 15),
                       stringsAsFactors = FALSE)

apply(data.in[76:125,], 1, function(x){
    custom.config <- umap.defaults
    custom.config$n_neighbors <- x[1]
    custom.config$min_dist <- x[2]
    custom.config$n_components <- x[3]
    
    fgwc_param <- c(kind = 'v', ncluster = 7, m = 1.1, distance = 'euclidean', 
                    order = 2, alpha = 0.5, a = 1, b = 1, max.iter = 500, 
                    error = 1e-5, randomN = 1)
    
    ## Run UMAP ----
    scaled <- "scaled"
    umap.custom <- umap::umap(inputUMAP.pca,
                              config = custom.config)
    
    ## Run FGWC ----
    umap.fgwc.in <- data.frame(umap.custom$layout[,])
    colnames(umap.fgwc.in) <- paste0("UMAP", 1:dim(umap.custom$layout)[2])
    
    umap.map.clst.FGWC <- naspaclust::fgwc(data = umap.fgwc.in, 
                                           pop = pop, 
                                           distmat = dist.Mat,
                                           algorithm = "classic",
                                           fgwc_param = fgwc_param)
    
    knn = custom.config$n_neighbors
    minDist = custom.config$min_dist
    lowD = custom.config$n_components
    col.No = length(unique(umap.map.clst.FGWC$cluster))
    colour.values <- get.colours(col.No)
    ### Set graph output prefix ----
    prefix.fgwc <- paste0("umap.custom.FGWC.")
    
    ### Plot graph on the tissue map ----
    umap.map <- umap.fgwc.in %>%
        mutate("geometry" = polygons$geom_pol) %>%
        mutate( "fgwc" = as.factor(umap.map.clst.FGWC$cluster))
    
    p <- ggplot() + 
        geom_sf(data = umap.map$geometry,
                aes(fill = umap.map$fgwc)) + 
        scale_fill_manual(values = colour.values) +
        labs(title = "FGWC clusters",
             caption = paste0("knn = ", custom.config$n_neighbors,
                              "\nminDist = ", custom.config$min_dist,
                              "\nlowDim = ", lowD,
                              "\nscaled = ", scaled),
             fill = "FGWC clusters") + 
        my_theme
    #print(p)
    ggsave(file.path(gwumapDir, paste0(prefix.fgwc, ".map.pca.", scaled, ".knn", knn, ".minD", minDist, ".lowDim", lowD, ".tiff")),
           plot = p,
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 150)
})


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# 10. Run GW-UMAP - weight gene expression ----
## make observations matrix g x s ----
obs <- counts[select,] %>% 
    .[,nb_names] %>%
    t() %>% 
    as.data.frame()

## make spatial weights matrix s x s ----
### Set parameters ----
dist.Mat <- gw.dist(dp.locat = st_coordinates(polygons$geom_cntd), p = 2)
bw = (range(dist.Mat)[2])+4
kernel = "bisquare"
### Calculate W distance matrix ----
w <- gw.weight(vdist = dist.Mat, 
               bw = bw,
               kernel = kernel,
               adaptive = FALSE)
### Plot a weighted distance example ----
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

## make the array s x g x s ----
### Set dimensions
loc.n <- nrow(obs)
var.n <- ncol(obs)
umap.n <- loc.n

obs.W <- array(data = 0, c(umap.n, loc.n, var.n))

## weight expression ----
### Select weights for focus point ----
focus = 1
w.focus <- w[focus,]
### Apply weights ----
sweep <- sweep(obs, 1, w.focus, "*")
sweep <- colSums(sweep)
sumW <- sum(w.focus)
sweep <- sweep(obs, 2, sweep/sumW)
sqrt.W <- sqrt(w.focus)
sweep <- sweep(sweep, 1, sqrt.W, "*")

## UMAP ----
### Set parameters ----
custom.config <- umap.defaults
custom.config$n_neighbors <- 7
custom.config$min_dist <- 0.1
custom.config$n_components <- 2
### Run UMAP ----
umap <- umap::umap(sweep,
                   config = custom.config)

umap.layout <- as.data.frame(umap$layout[,])

# ------------------------------------------- #
# ------------------------------------------- #
umap.knn <- umap$knn$indexes %>%
    as.numeric() %>%
    matrix(nrow = 1185, ncol = 7)

umap.dist <- umap$knn$distances %>%
    as.numeric() %>%
    matrix(nrow = 1185, ncol = 7)

dist.sweep <- dist(sweep) %>% 
    as.matrix()

n <- nrow(dist.sweep)
k <- 7
dist.sweep.knn <- matrix(0, ncol = k, nrow = n)
dist.sweep.dist <- knn.mat
for(i in 1:n){
    dist.sweep.knn[i,] = order(dist.sweep[i,])[1:k]
    dist.sweep.dist[i,] = dist.sweep[i,dist.sweep.knn[i,]]
}

dist.sweep.knn
dist.sweep.dist

rm(knn.mat, n, k, knd.mat)

identical(umap.knn, dist.sweep.knn)
identical(umap.dist, dist.sweep.dist)
# ------------------------------------------- #
# ------------------------------------------- #
fgwc_param <- c(kind = 'u', ncluster = 7, m = 1.1, distance = 'euclidean', 
                order = 2, alpha = 0.5, a = 1, b = 1, max.iter = 500, 
                error = 1e-5, randomN = 1)
opt_param <- c(vi.dist='uniform',npar=10,par.no=2,par.dist='euclidean',par.order=2,pso=TRUE,
               same=10,type='sim.annealing',ei.distr='normal',vmax=0.7,wmax=0.9,wmin=0.4,
               chaos=4,x0='F',map=0.7,ind=1,skew=0,sca=1)

umap.fgwc.in <- umap.layout
colnames(umap.fgwc.in) <- paste0("UMAP", 1:dim(umap.fgwc.in)[2])

umap.map.clst.FGWC <- naspaclust::fgwc(data = umap.fgwc.in, 
                                       pop = pop, 
                                       distmat = dist.Mat,
                                       algorithm = "abc",
                                       fgwc_param = fgwc_param,
                                       opt_param = opt_param)

umap.map <- umap.fgwc.in %>%
    mutate("geometry" = polygons$geom_pol) %>%
    mutate( "fgwc" = as.factor(umap.map.clst.FGWC$cluster))

ggplot() + 
    geom_sf(data = umap.map$geometry,
            aes(fill = umap.map$fgwc)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "FGWC clusters",
         caption = paste0("knn = ", custom.config$n_neighbors,
                          "\nminDist = ", custom.config$min_dist,
                          "\nlowDim = ", lowD,
                          "\nscaled = ", scaled),
         fill = "FGWC clusters") + 
    my_theme
# ------------------------------------------- #
# ------------------------------------------- #
umap.hdbsc.in <- umap.layout
colnames(umap.hdbsc.in) <- paste0("UMAP", 1:dim(umap.hdbsc.in)[2])
umap.map.clst.HDBSC <- hdbscan(umap.hdbsc.in, 
                               minPts = 12,
                               gen_hdbscan_tree = TRUE,
                               verbose = TRUE)

plot(umap.map.clst.HDBSC)

minPts = umap.map.clst.HDBSC$hc$call$minPts
knn = custom.config$n_neighbors
minDist = custom.config$min_dist
lowD = custom.config$n_components
col.No = length(unique(umap.map.clst.HDBSC$cluster))
colour.values <- get.colours(col.No)
### Set graph output prefix ----
prefix.hdbscan <- paste0("umap.custom.HDBSCAN.")

### Plot graph on the tissue map ----
umap.map <- umap.hdbsc.in %>%
    mutate("geometry" = polygons$geom_pol) %>%
    mutate( "hdbscan" = as.factor(umap.map.clst.HDBSC$cluster))

ggplot() + 
    geom_sf(data = umap.map$geometry,
            aes(fill = umap.map$hdbscan)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "HDBSCAN clusters",
         caption = paste0("knn = ", custom.config$n_neighbors,
                          "\nminDist = ", custom.config$min_dist,
                          "\nminPts = ", minPts,
                          "\nlowDim = ", lowD,
                          "\nscaled = ", scaled),
         fill = "HDBSCAN clusters") + 
    my_theme
