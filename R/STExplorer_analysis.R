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
library(cols4all)

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
    my_theme

ggsave(file.path(graphDir, "voronoi_tessellation_no.colour.tiff"),
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
    geom_sf(data = nb_sf, colour = "black") + 
    geom_point(data = polygons, aes(x = pixel_x, y = pixel_y, colour = factor(nb_count))) + 
    # Add titles and visually format the plot:
    scale_color_manual(values = c("#34568B", "#FF6F61", "#88B04B",
                                  "#FDAC53", "#F7CAC9", "#6B5B95")) +
    labs(title = paste("Contiguity neighbours"),
         subtitle = "",
         colour = "Neighbours\ncount") + 
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme

ggsave(file.path(graphDir, "voronoi_tessellation_on-tissue.tiff"),
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

# Get spot names
nb_names <- polygons$Barcode

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
# Generate a df with combinations of parameters
data.in <- param.combo(var.no = c(500, 750, 1000),
                       k = c(20, 50, 100),
                       kernel = c("gaussian", "exponential"))

# Initialise an empty list
pca_gw.list <- list()

# Run multiple GWPCAs and output them into the list
pca_gw.list <- gwpca.combo(data.in)

# Get the parameters in a table
pca_params <- cbind(sapply(pca_gw.list, get.params)) %>%
    t() %>%
    as.data.frame() %>%
    mutate(minutes = round(as.numeric(minutes), 2))

## Plot global PCA results ----
# Prepare data for scree plot
pvar <- pca_gw.list$pca_gw.500.20.gau$pca$sdev^2/sum(pca_gw.list$pca_gw.500.20.gau$pca$sdev^2)
pvar <- data.frame(var = pvar,
                   PCs = sprintf("PC%02d", seq(1, 500)))

# Plot scree plot
ggplot(pvar[1:10,], aes(x = PCs, y = var, group = 1)) + 
    geom_point(size=3)+
    geom_line() +
    xlab("Principal Component") +
    ylab("% Variance Explained") +
    ggtitle("Scree Plot") +
    ylim(0, 1) + 
    my_theme

ggsave(file.path(graphDir, "gwpca_globalPCA_screeplot.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Find leading items at each location
lead.item <- gwpca.Leading.G.single(pca_gw.list$pca_gw.500.20.gau, 
                                    pc.no = 1, 
                                    sf.geom = polygons$geom_pol,
                                    gene.names = TRUE,
                                    biomart = biomart.mouse.98,
                                    check.names = FALSE) 

pc.No = 1
col.No = length(unique(lead.item[,1]))

ggplot() + 
    geom_sf(data = lead.item$geometry,
            aes(fill = lead.item[,1])) +
    scale_fill_manual(values = c4a("wright25", col.No)) + 
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    labs(title = paste0("Leading Genes on PC", pc.No),
         fill = "Leading Genes") +
    my_theme + 
    theme(legend.position="right")

ggsave(file.path(graphDir, 
                 paste0("gwpca_.500.20.gau_leadingGene_PC", pc.No, ".tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Find top leading genes in each location 
dt_top.Gs <- gwpca.topLeading.Gs(gwpca = pca_gw.list$pca_gw.500.20.gau,
                                 pc = 1,
                                 genes.n = 3,
                                 sf.geom = polygons$geom_pol,
                                 method = "membership",
                                 gene.names = TRUE,
                                 biomart = biomart.mouse.98,
                                 check.names = FALSE)


pc = 1
genes = 3
method = "membership"
groups = count(unique(as.data.frame(dt_top.Gs$Top_lead_Gs)))

#map the leading genes
ggplot() +
    geom_sf(data = dt_top.Gs$geometry, 
            aes(fill = as.factor(dt_top.Gs$Top_lead_Gs))) +
    geom_point(size = 3) + 
    #scale_fill_manual(values = c4a("wright25", 18)) +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    labs(title = paste0("Leading Genes on PC", pc),
         subtitle = paste0("Top ", genes, " Genes"),
         caption = paste0("Grouping method: ",
                          method,
                          "\nNumber of groups: ",
                          groups),
         fill = "Top Leading Genes\nGroups") +
    my_theme +
    theme(legend.position="right")

ggsave(file.path(graphDir, "gwpca_.500.20.gau_leadingGenes3_PC1-membership.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Calculate the PTV for multiple Components
props <- gwpca.prop.var(gwpca.obj = pca_gw.list$pca_gw.500.20.gau,
                        n.comp = c(5, 10, 20, 30, 40, 50))

# Plot them all together
gwpca.plot.prop.vars.multi(select(props, -c(pixel_x, pixel_y)), theme = my_theme)

ggsave(file.path(graphDir, "gwpca_.500.20.gau_PTV_all.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Plot PTVs one by one and in a panel
## One by one:
for (n in column_names) {
    gwpca.plot.prop.vars.single(data = select(props, n))
    
    name = sub("Comps_", "", n)
    
    ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_PTV_", name,".pdf")),
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 400)
}

## In a panel
plot_list <- lapply(colnames(select(props, -c(pixel_x, pixel_y))),
                    gwpca.plot.prop.vars.single, 
                    data = select(props, -c(pixel_x, pixel_y)), ylab = NULL)

plot_list <- setNames(plot_list, colnames(select(props, -c(pixel_x, pixel_y))))

egg::ggarrange(plots = plot_list, nrow = 2, ncol = 3)

ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_PTV_panel.pdf")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Map the PTV for a specific selection of components
for (i in c(5, 10, 20, 30, 40, 50)) {
    comps <- sprintf("Comps_%02d", i)
    ptv.map <- dplyr::select(props, all_of(c(comps, "geometry")))
    
    ggplot() + 
        geom_sf(data = ptv.map$geometry,
                aes(fill = ptv.map[,1])) +
        scale_fill_viridis_c(option = "magma", limits = c(0, 100)) +
        xlab("X coordinates (pixels)") +
        ylab("Y coordinates (pixels)") +
        labs(title = "Percantage of Total Variation\n(PTV)",
             fill = paste0("PTV of ", i, "\n components")) +
        my_theme
}

# calculate the discrepancies
data.mat <- as.matrix(inputPCAgw@data)
discrepancy <- gwpca.cv.contrib(data.mat, coordinates(inputPCAgw), 
                                bw =6*spot_diameter(spatialDir), 
                                adaptive = TRUE, dMat = dist.Mat)

# plot the discrepancies in a box plot
discrepancy_df <- data.frame(disc = discrepancy)

ggplot(pivot_longer(discrepancy_df, col = "disc"),
       aes(x = name, y = value)) + 
    geom_boxplot(fill = "#D1E5F0", colour = "#2166AC", 
                 outlier.colour = "red", outlier.size = 2) + 
    geom_jitter(col = "#EF8A62", size = 2, width = 0.3, alpha = 0.8) +
    coord_flip() +
    ggtitle("Local PC Discrepancy") +
    xlab(NULL) + 
    my_theme

ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_discreps.box.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# map the discrepancies
dt <- inputPCAgw@data %>%
    mutate(disc = discrepancy,
           geometry = polygons$geom_pol)

disc.map <- dplyr::select(dt, all_of(c("disc", "geometry")))

ggplot() + 
    geom_sf(data = disc.map$geometry, 
            aes(fill = disc.map$disc)) + 
    scale_fill_viridis_c(option = "inferno") +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    labs(title = "Local PC Discrepancy",
         fill = "Discrepancy\nscore") +
    my_theme

ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_discreps.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# create the input data table for the outlier plot
inputPCAgw.outlier <- vst_df[,select] %>% # select top 500 variable genes
    as.data.frame() %>%                      # make it a df
    .[nb_names,]                             # order rows

# create the biomart for mouse
biomart.mouse <- create_biomart("mouse")

# plot the heatmap to visualise the genes that make this location an outlier
tiff(file.path(graphDir, paste0("gwpca_.500.20.gau_discrep.max.outlierPlot.tiff")),
     width = grDevices::dev.size(units = "in")[1],
     height = grDevices::dev.size(units = "in")[2],
     units = "in",
     res = 400)

gwpca.plot.outlier(inputPCAgw.outlier,
                   bw = 3*spot_diameter(spatialDir),
                   focus = which(discrepancy_df == max(discrepancy_df)),
                   dMat = dist.Mat,
                   show.vars = "top",
                   mean.diff = 1,
                   gene.names = TRUE, 
                   biomart = biomart.mouse.98,
                   show.data = FALSE,
                   check.names = FALSE,
                   scale = "row",
                   cutree_cols = 5,
                   cutree_rows = 6,
                   color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(1000)))

dev.off()

## Prepare for Fuzzy Geographically Weighted Clustering (FGWC) ----
# Calculate the weighted distance matrix
dist.Mat<- gw.dist(dp.locat = coordinates(inputPCAgw))

# Generate a population matrix
pop <- as.matrix(rep(1, nrow(vst_df)))

# Select only the top variable genes to drive the clustering
inputFGWC <- vst_df[select] %>% # select 500 most variable genes
    .[nb_names,] # Re-order rows to match the polygon object row order

# Set FGWC parameters
fgwc_param <- c(kind = 'v', ncluster = 7, m = 1.1, distance = 'euclidean', 
                order = 2, alpha = 0.5, a = 1, b = 1, max.iter = 500, 
                error = 1e-5, randomN = 1)

## Run FGWC ----
fgwc <- naspaclust::fgwc(data = inputFGWC, 
                         pop = pop, 
                         distmat = dist.Mat,
                         algorithm = "classic",
                         fgwc_param = fgwc_param)

fgwc_clusters <- data.frame(geometry = polygons$geom_pol,
                            cluster = fgwc$cluster)

ggplot() + 
    geom_sf(data = fgwc_clusters$geometry, 
            aes(fill = as.factor(fgwc_clusters$cluster)),
            colour = "grey30", 
            show.legend = TRUE) +
    scale_fill_brewer(type = "qual") + 
    labs(title = "Fuzzy GW Clustering (FGWC)",
         subtitle = "nclust = 7",
         fill = "Cluster") + 
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme

ggsave(file.path(graphDir, "fgwc_nclust-8.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

#---------------------TEST STUF...------------------------------#
#---------------------------------------------------------------#

