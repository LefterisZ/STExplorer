# This script will deal with extracting information needed for neighbour
#   identification and weighted matrix calculation from the Spaniel sce
#   object.

# packages needed
library(spdep)
library(sf)
library(jsonlite)
library(ggplot2)
library(naspaclust)
library(rlist)

# non-package functions needed
source("./R/spot_diameter.R")
source("./R/sf_coord_as_df.R")
source("./R/sfc_coord_as_df.R")
source("./R/readSpacerangerMD.R")
source("./R/readSpacerangerD.R")

## set the file paths to spaceranger's spatial and gene expression folders ----
sampleDir <- "Olfactory_Bulb/Olfactory_Bulb_A1_Results"
spatialDir <- file.path(inputDir, sampleDir, "spatial")
countsDir <- file.path(inputDir, sampleDir, "filtered_feature_bc_matrix")

## Import the dataset ----
inputMD <- readSpacerangerMD(spatialDir, res = "low") #read-in MetaData
inputD <- readSpacerangerD(countsDir) #read-in gene expression Data

## Select spots in both bins (Sections) 0 and 1 ----
spot_position <- inputMD %>% 
    select(c("Barcode", "pixel_x", "pixel_y", "Section"))

## Convert spots to centroids ----
centroids <- spot_position %>% 
  st_as_sf(coords = c("pixel_x", "pixel_y"), 
           remove = FALSE)

## Combine the points into a multipoint geometry: ----
cntd_union <- st_union(centroids)
head(cntd_union)

## Use the union of points to generate a voronoi object ----
voronoi <- st_voronoi(cntd_union, bOnlyEdges = TRUE)
head(voronoi)

## Create an enveloped voronoi tessellation around the tissue ----
voronoi_env <- st_intersection(st_cast(voronoi), st_convex_hull(cntd_union))
head(voronoi_env)

## plot the voronoi tessellation ----
ggplot() +
    geom_sf(data = voronoi_env, colour = "black", fill = "white") + 
    geom_sf(data = centroids, aes(colour = as.factor(Section)), size = 1.2) + 
    # Add titles and visually format the plot:
    labs(title = paste("Voronoi tessellation"),
         subtitle = ,
         colour = "Spot\nPosition") + 
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") +
    scale_colour_manual(values = c("0" = "#EEA47FFF", "1" = "#00539CFF"),
                        labels = c("off-tissue", "on-tissue")) +
    theme(axis.title = element_text(size = rel(2)),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = rel(2)),
          axis.ticks.length.y.left = unit(.15, "cm"),
          axis.ticks.length.x.bottom = unit(.15, "cm"),
          plot.title = element_text(colour = "black", size = rel(2.5), hjust = 0.5),
          plot.subtitle = element_text(colour = "black", size = rel(1.9), hjust = 0.5),
          plot.margin = unit(c(0.5,1,0.5,1), "cm"),
          legend.title = element_text(colour = "black", size = rel(1.5)),
          legend.text = element_text(colour = "black", size = rel(1.2)),
          panel.background = element_rect(fill = "white"))


ggsave(file.path(graphDir, "voronoi_tessellation.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Generate the POLYGONS from the MULTILINESTRING for bin_1 only and attach the ----  
##  barcode names
polygons <- st_polygonize(voronoi_env) %>% # polygonise the tessellation
    st_cast() %>% # convert GEOMETRYCOLLECTION to multiple POLYGONS
    st_sf() %>%  # convert sfc object to sf for st_join afterwards
    st_join(., 
            centroids[centroids$Section == 1,],
            join = st_contains,
            left = FALSE) %>% # Join the centroids with the POLYGONS
    mutate(Barcode_rn = Barcode) %>% # duplicate the barcode column
    column_to_rownames("Barcode_rn") %>% # move duplicate column to row names
    st_sf() # convert back to sf (mutate makes it a df)

## Create contiguity neighbours ----
neighbours <- poly2nb(polygons, snap = 0)
names(neighbours) = attr(neighbours, "region.id") # add names to the sub-lists

## Add number of neighbours for each polygon back to the polygons object ----
polygons$nb_count <- card(neighbours)

## Add the neighbour (nb) IDs as a nested df in the polygons object ----
nb_IDs <- neighbours %>%
    nb2lines(., coords = polygons$geometry) %>% #get nb connecting lines 
    as("sf") %>% #convert to sf
    st_drop_geometry() %>% #drop geometry column
    select(i_ID, j_ID) %>% #select only nb ID columns
    rename(nb_IDs = j_ID) %>% #rename the neighbours ID column
    group_by(i_ID) %>% #group by spot
    nest() #nest the groupings

polygons <- right_join(polygons, nb_IDs, by = c("Barcode" = "i_ID")) %>%
    rename(nb_IDs = data, geom_pol = geometry)

## Update the polygon object to keep the centroid geometries as well ----
polygons <- left_join(as.data.frame(polygons), as.data.frame(centroids), 
                      by = c("Barcode" = "Barcode"), suffix = c("", ".y")) %>%
    select(!ends_with(".y")) %>% 
    rename(geom_cntd = geometry) %>% 
    st_sf(sf_column_name = "geom_pol")


## Get a neighbours object for ggplot2 plotting ----
nb_sf <- as(nb2lines(neighbours, coords = polygons$geom_cntd), "sf")

## Plot neighbours graph ----
ggplot() +
    geom_sf(data = polygons$geom_pol, colour = "grey30", fill = "white") +
    #geom_sf(data = nb_sf, colour = "black") + 
    #geom_point(data = polygons, aes(x = pixel_x, y = pixel_y, colour = factor(nb_count))) + 
    geom_circle(data = data_circle, aes(x0 = x0, y0 = y0, r = r)) +
    # Add titles and visually format the plot:
    scale_color_manual(values = c("#34568B", "#FF6F61", "#88B04B",
                                  "#FDAC53", "#F7CAC9", "#6B5B95")) +
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

ggsave(file.path(graphDir, "voronoi_polygons_only.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

## Calculate neighbour weights with a distance decay function ----
neighbours_wght <- nb2listwdist(neighbours, polygons$geom_cntd,
                                type = "idw", style = "raw", alpha = 1)

## Import gene counts ----
inputD <- readSpacerangerD(countsDir)

## Prepare for gene expression normalisation using DESeq2 ----
spotName <- colnames(inputD)
spotTable <- data.frame(spotName = spotName)

#filter out genes with less than 10 counts total in all spots
inputD_filt <- inputD[!rowSums(inputD) < 10,]

dds <- DESeqDataSetFromMatrix(countData = inputD_filt,
                              colData = spotTable,
                              design = ~spotName)

dds = estimateSizeFactors(dds) # Estimate size factors

counts = counts(dds, normalized = TRUE) # export normalised counts

## Prepare for Geographically Weighted PCA (GWPCA) ----
# Transform counts to vst
vst <- varianceStabilizingTransformation(dds)
vst_df <- as.data.frame(t(assay(vst))) # transpose and transform to df

# Get the coordinates
coords <- polygons[, c("Barcode", "pixel_x", "pixel_y")] %>%
    st_drop_geometry() %>% 
    column_to_rownames(var = "Barcode")

# Get the data into a SpatialPointsDataFrame object
inputPCAgw <- SpatialPointsDataFrame(coords, vst_df, match.ID = TRUE)

# Identify the most variable genes equal to the number of spots.
# gwpca uses princomp to run the PCAs and this does not accept the number of
# variables (genes) being more than the number of samples (spots).
row_vars <- rowVars(assay(vst))
select <- order(row_vars, decreasing = TRUE)[seq_len(500)]
inputPCAgw <- inputPCAgw[select]
vars <- colnames(inputPCAgw@data)
bw <- 6*spot_diameter(spatialDir)
k <- 20

## Run GWPCA ----
pca_gw <- gwpca(inputPCAgw, 
                vars = vars, 
                bw = bw,
                k = k,
                kernel = "gaussian")

#### RUN MULTIPLE GWPCAs ####
#function to generate a parameters combo df
param.combo <- function(...){
    data <- expand.grid(..., stringsAsFactors = FALSE)
    data$kernl <- substring(data$kernel, 1, 3)
    data$obj <- paste0("pca_gw.", data$var, ".", data$k, ".", data$kernl)
    data$kernl <-NULL
    print(head(data, 3))
    
    return(data)
}

# function to run gwpca
gwpca.param.combo <- function(var.no, k, kernel){
    inputPCAgw <- SpatialPointsDataFrame(coords, vst_df, match.ID = TRUE)
    print(var.no)
    select <- order(row_vars, decreasing = TRUE)[seq_len(var.no)]
    inputPCAgw <- inputPCAgw[select]
    vars <- colnames(inputPCAgw@data)
    bw <- 6*spot_diameter(spatialDir)
    k <- k
    obj <- gwpca(inputPCAgw, 
                 vars = vars, 
                 bw = bw,
                 k = k,
                 kernel = kernel)
    
    return(obj)
}

# function that wraps gwpca.param.combo and outputs objects in a list
gwpca.combo <- function(var.no, k, kernel, list) {
    list <- list.append(list, gwpca.param.combo(var.no, k, kernel))
    return(list)
}

# function to extract parameters
get.params <- function(gwpca.out, obj){
    obj <- data.frame("vars" = length(gwpca.out$GW.arguments$vars),
                      "spots" = gwpca.out$GW.arguments$dp.n,
                      "PCs" = gwpca.out$GW.arguments$k,
                      "kernel" = gwpca.out$GW.arguments$kernel,
                      "minutes" = gwpca.out$timings$stop - gwpca.out$timings$start)
    
    return(obj)
}


data.in <- param.combo(var.no = c(500, 750, 1000),
                       k = c(20, 50, 100),
                       kernel = c("gaussian", "exponential"))

pca_gw.list <- list()

pca_gw.list <- with(data.in,
                    mapply(gwpca.combo, 
                           var.no, k, kernel, 
                           MoreArgs = list(list = pca_gw.list)))

pca_gw.list <- setNames(pca_gw.list, data.in$obj)


## Prepare for Fuzzy Geographically Weighted Clustering (FGWC) ----
# Calculate the weighted distance matrix
dist.Mat<- gw.dist(dp.locat = coordinates(inputPCAgw), 
                   rp.locat = coordinates(inputPCAgw))

# Generate a population matrix
pop <- as.data.frame(rep(1, nrow(vst_df)))

# Select only the top variable genes to drive the clustering
inputFGWC <- vst_df[select]

# Set FGWC parameters
fgwc_param <- c(kind = 'u', ncluster = 10, m = 1.2, distance = 'euclidean', 
                order = 2, alpha = 0.5, a = 1, b = 1, max.iter = 500, 
                error = 1e-5, randomN = 1)

## Run FGWC ----
fgwc <- naspaclust::fgwc(data = inputFGWC, 
                         pop = pop, 
                         distmat = dist.Mat,
                         algorithm = "classic",
                         fgwc_param = fgwc_param)

#---------------------TEST STUF...------------------------------#
#---------------------------------------------------------------#

