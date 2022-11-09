install.packages("umap")
install.packages("graphlayouts")
install.packages("ggraph")
install.packages("snahelper")
library(umap)
library(graphlayouts)
library(igraph)
library(ggraph)

# library(plotly)

rm(E,xy)

# 1. Create a gwumap graphics output directory ----
gwumapDir <- file.path(graphDir, "gwumap")

# 2. Create input data from 500 most variable genes ----
inputUMAP <- vst_df[, select] %>%
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

ggraph::ggraph(g = umap.map.3.ig,
               layout = "stress") + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = grp)) + 
    scale_colour_manual(values = c4a("wright25", 7))

### Cluster with Louvain ----
umap.map.3.clst <- cluster_louvain(umap.map.3.ig)
membership(umap.map.3.clst)
plot(umap.map.3.clst, umap.map.3.ig)
V(umap.map.3.ig)$louv <- as.factor(membership(umap.map.3.clst))

ggraph::ggraph(g = umap.map.3.ig,
               layout = "stress") + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = c4a("wright25", 9))

ggraph::ggraph(g = umap.map.3.ig,
               layout = "manual",
               x = umap.custom.3$layout[,1],
               y = umap.custom.3$layout[,2]) + 
    geom_edge_link(width = 0.1, edge_colour = "grey66") + 
    geom_node_point(aes(colour = louv)) + 
    scale_colour_manual(values = c4a("wright25", 9))

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
kernel = "gaussian"

polygons.wdist <- gw.weight(vdist = dist.Mat, 
                            bw = bw,
                            kernel = kernel,
                            adaptive = FALSE)
range(polygons.wdist)
range(dist.Mat)

### Calculate weighted UMAP distances
umap.map.3.wdist <- umap.map.3.dist * dist.Mat

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
         subtitle = "kernel = bisquared") +
    my_theme

### Reconstruct tissue map from simple point distances ----
cmdscale_polygons <- cmdscale(dist.Mat, k = 2, add = TRUE) %>%
    .[["points"]] %>%
    as.data.frame() %>%
    mutate(group = umap.map.3$group) %>%
    rename("X" = "V1", "Y" = "V2")
ggplot() +
    geom_point(data = cmdscale_polygons,
               aes(x = X, y = Y, colour = group), size = 4.5) +
    scale_color_manual(values = c4a("wright25", 7)) +
    labs(title = "Reconstructed map from point distances") +
    my_theme

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
