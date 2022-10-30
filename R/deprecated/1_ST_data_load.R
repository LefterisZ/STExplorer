######################################################-
# This is a script to load,  cluster and save        #
#   Spatial Transcriptomics data                     #
#                                                    #
#' @author Eleftherios (Lefteris) Zormpas            #
#                                                    #
# This script is a replicate/ adaptation from an     #
#   original script written and provided by          #
#   Dr Rachel Queen Newcastle University             #
#   ---> give her the credit!                        #
######################################################-

#----------------------------------------------------#
# START OF SCRIPT ----
#----------------------------------------------------#


#----------------------------------------------------#
# 1. SET FILE PATHS ----
#----------------------------------------------------#
## Input data folder ----
dataDir <- file.path(inputDir, sampleName)

## Output data folder ----
sce_output <- paste0(outputDir, sampleName, "_sce.rds")
batch_corrected_seurat_output <- paste0(outputDir, sampleName, "_seurat.rds")


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
  sce_list[[i]] <- Spaniel::createVisiumSCE(file.path(dataDir, filefolders[i], "Olfactory_Bulb_A1_Results"))
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

## 4.2 combine replicates to one object
seuratObj_combined <- merge(seurat_list[[1]], seurat_list[-1])
seuratObj_batch_corrected <- seuratObj_combined %>%
  ## 4.3 normalise for cell to cell differences
  Seurat::NormalizeData(normalization.method = "LogNormalize",
                        scale.factor = 10000,
                        verbose = FALSE) %>%
  ## 4.4 find genes to use for clustering
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ## 4.5 scale genes so level of expression is the similar for all genes
  ScaleData(verbose = FALSE,
            vars.to.regress = c("detected",
                                "total"))


seuratObj_batch_corrected <- seuratObj_batch_corrected %>%
  ## 4.6 dimension reduction --> linear
  RunPCA(
    features = VariableFeatures(seuratObj_batch_corrected),
    npcs = 20,
    verbose = FALSE
  ) %>%
  ## 4.7 batch correction (use group.by.var = "id" --> from meta.data -is the sample ID-)
  RunHarmony("id", plot_convergence = TRUE, assay.use = "originalexp")


#----------------------------------------------------#
# 5. SEURAT CLUSTERING ----
#----------------------------------------------------#
seuratObj_batch_corrected <- seuratObj_batch_corrected %>%  
  ## 5.1 Find Nearest Neighbours
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  ## 5.2 Cluster the data
  FindClusters(resolution = c(0.5)) %>%
  # run non-linear dimension reduction for later visualisation
  RunUMAP(reduction = "harmony", dims = 1:10, assay = "originalexp")

## 5.4 save RDS
saveRDS(seuratObj_batch_corrected, batch_corrected_seurat_output)

## 5.5 add clusters to sce object
ids <- sampleInfo$id
for (i in 1:4) {
  id <- ids[i]
  sce_list[[i]]$cluster <-
    seuratObj_batch_corrected[, seuratObj_batch_corrected$id == id]$originalexp_snn_res.0.5
}

## 5.6 save RDS again with added clusters
saveRDS(sce_list, sce_output)

## 5.7 plot the clusters
seurat_clust <- as.data.frame(colData(sce_list$Olfactory_Bulb))
ggplot(seurat_clust) + 
    geom_point(aes(x = pixel_x, y = pixel_y, colour = cluster), size = 3) + 
    scale_colour_brewer(type = "qual") +
    scale_y_reverse() +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    ggtitle("Seurat Clustering") +
    my_theme

ggsave(file.path(graphDir, "seurat_clustering_mob.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

#----------------------------------------------------#
# END OF SCRIPT ----
#----------------------------------------------------#