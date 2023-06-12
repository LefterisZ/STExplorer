
#----------------------------------------------------#
# 1. SET FILE PATHS ----
#----------------------------------------------------#
## Input data folder ----
sampleName  <- list.files(path = inputDir)
dataDir <- file.path(inputDir, sampleName)

#----------------------------------------------------#
# 2. CREATE AN INFO TABLE ----
#----------------------------------------------------#
## Create a sampleInfo table ----
filefolders <- list.files(dataDir)
id <- filefolders %>% gsub("results_", "", .) %>% gsub("_Results", "", .)
section <- id %>% gsub(".*_", "", .)
sampleInfo <- data.frame(
    fileFolders = filefolders,
    id = id,
    tissue = sampleName,
    section = section
)

## Double-check the sampleInfo table
sampleInfo

#----------------------------------------------------#
# 3. GENERATE AN SCE LIST ----
#----------------------------------------------------#
## Create an sce list (SCE: SingleCellExperiment Object)
sce_list <- list()
for (i in 1:length(filefolders)) {
    ## create a VisiumSCE object - defaults to low res image!!
    sce_list[[i]] <- Spaniel::createVisiumSCE(file.path(dataDir, filefolders[i]))
    ## add QC metrics for low quality spots identification
    sce_list[[i]] <- addPerCellQC(sce_list[[i]], subsets = list(Mito = 1:10))
    ## add info from sampleInfo to colData
    sce_list[[i]]$id <- sampleInfo$id[i]
    sce_list[[i]]$tissue <- sampleInfo$tissue[i]
    sce_list[[i]]$section <- sampleInfo$section[i]
}

## add the names from SCE list into the sampleInfo
names(sce_list) <- sampleInfo$fileFolders

#----------------------------------------------------#
# 4. DATA PRE-PROCESSING ----
#----------------------------------------------------#
## 4.1 create seurat list from the sce list
seurat_list <- lapply(sce_list, as.Seurat, data = NULL)

## 4.2a combine replicates to one object
seuratObj_combined <- merge(seurat_list[[1]], seurat_list[-1])

## 4.2.b if no multiple replicates exist:
seuratObj_combined <- seurat_list[[1]]

seuratObj_batch_corrected2 <- seuratObj_combined %>%
    ## 4.3 normalise for cell to cell differences
    Seurat::NormalizeData(normalization.method = "LogNormalize",
                          scale.factor = 10000,
                          verbose = TRUE) %>%
    ## 4.4 find genes to use for clustering
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ## 4.5 scale genes so level of expression is the similar for all genes
    ScaleData(verbose = TRUE,
              vars.to.regress = c("detected",
                                  "total"))


seuratObj_batch_corrected <- seuratObj_batch_corrected2 %>%
    # 4.6 dimension reduction --> linear
    RunPCA(features = VariableFeatures(seuratObj_batch_corrected2),
           npcs = 20,
           verbose = TRUE) %>%
    ## 4.7 batch correction (use group.by.var = "id" --> from meta.data -is the sample ID-)
    # run the below only if you have multiple slices combined. If not do not run
    RunHarmony("id", plot_convergence = TRUE, assay.use = "originalexp")

saveRDS(seuratObj_batch_corrected_norm, file = paste0(outputDir, "DLPFC_151673_counts_Seurat_Norm.rds"))
saveRDS(seuratObj_batch_corrected2, file = paste0(outputDir, "DLPFC_151673_counts_Seurat_Norm_VST_Var2000_Scaled.rds"))
saveRDS(seuratObj_combined, file = paste0(outputDir, "DLPFC_151673_SeuratObject.rds"))

#----------------------------------------------------#
# 5. SEURAT CLUSTERING ----
#----------------------------------------------------#
# If you run harmony above then use reduction = "harmony". If not then use
# reduction = "pca"
seuratObj_batch_corrected <- seuratObj_batch_corrected %>%  
    ## 5.1 Find Nearest Neighbours
    FindNeighbors(reduction = "pca", dims = 1:10) %>%
    ## 5.2 Cluster the data
    FindClusters(resolution = c(0.8)) %>%
    # run non-linear dimension reduction for later visualisation
    RunUMAP(reduction = "pca", dims = 1:10, assay = "originalexp")

## 5.4 save RDS
#saveRDS(seuratObj_batch_corrected, batch_corrected_seurat_output)

## 5.5 add clusters to sce object
ids <- sampleInfo$id
for (i in 1:4) {
    id <- ids[i]
    sce_list[[i]]$cluster <-
        seuratObj_batch_corrected[, seuratObj_batch_corrected$id == id]$originalexp_snn_res.0.5
}

## 5.6 save RDS again with added clusters
#saveRDS(sce_list, sce_output)

## 5.7 plot the clusters
# If transferred to the sce_list above run this:
seurat_clust <- as.data.frame(colData(sce_list$Olfactory_Bulb))
# If not, then run this:
seurat_clust <- as.data.frame(seuratObj_batch_corrected@meta.data)
seurat_clust$Barcode <- gsub("-1", ".1", seurat_clust$Barcode)

map <- seurat_clust %>% 
    left_join(polygons[,c("Barcode", "geom_pol")])

col.No = length(unique(map$originalexp_snn_res.0.75))
colour.values <- get.colours(col.No)

ggplot(seurat_clust) + 
    geom_sf(data = map$geom_pol,
            aes(fill = map$originalexp_snn_res.0.75)) + 
    scale_fill_manual(values = colour.values) +
    labs(title = "Seurat clustering",
         fill = "Clusters") + 
    my_theme

ggsave(file.path(graphDir, "seurat_clustering_mob_res075.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

seurat_clust_out <- dplyr::select(seurat_clust, c(Barcode, cluster))
write.csv(seurat_clust_out, file = "clusters_seurat.csv", row.names = FALSE)


## 5.8 get the scaled/normalised count data in a df
inputD_S <- as.data.frame(seuratObj_batch_corrected@assays$originalexp@scale.data)
colnames(inputD_S) <- seurat_clust$Barcode

## 5.9 get the normalised counts only in a df
inputD_S.norm <- seuratObj_batch_corrected@assays$originalexp@data %>% 
    as.matrix() %>% 
    as.data.frame()
colnames(inputD_S.norm) <- seurat_clust$Barcode

