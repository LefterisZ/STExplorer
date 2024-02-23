sampleDir_p <- c("./data/test_data/Visium_Human_Prostate/Patient1/H1_4",
               "./data/test_data/Visium_Human_Prostate/Patient1/H1_5",
               "./data/test_data/Visium_Human_Prostate/Patient1/V1_2")

sampleNames_p <- c("H1_4", "H1_5", "V1_2")

names(sampleDir_p) <- sampleNames_p

sampleNames_p

## Create the MSFE object
msfe_p <- MetaSpatialFeatureExperiment()

for (i in seq_along(sampleNames_p)) {
  message("Adding sample: ", sampleNames_p[i])
  msfe_p <- addSFE(msfe_p,
                 read10xVisiumSFE(samples = sampleDir_p[i],
                                  sample_id = sampleNames_p[i],
                                  type = "HDF5",
                                  data = "filtered",
                                  images = "lowres",
                                  style = "W",
                                  zero.policy = TRUE))
}

for (id in sampleNames_p) {
  message("Working on sample: ", id)
  ## Add location-related statistics
  msfe_p <- addPerLocQC(msfe_p,
                      sample_id = id,
                      gTruth = gTruth_list[[id]],
                      assay = "counts",
                      MARGIN = 2,
                      subsets = list(mito = is_mito[[id]]))
  message("\tAdded location-related statistics")

  ## Add geometries
  msfe_p <- addGeometries(msfe_p,
                        samples = sampleDir_p[id],
                        sample_id = id,
                        res = "fullres")
  message("\tAdded geometries")

  ## Add gene/feature-related statistics
  msfe_p <- addPerGeneQC(msfe_p,
                       sample_id = id,
                       assay = "counts",
                       version = NULL,
                       mirror = NULL)
  message("\tAdded gene/feature-related statistics")
}

capt_area_p <- list()
for (i in sampleNames_p) {
  data_p <- read.csv(file.path(sampleDir_p[i], "outs/spatial",
                                    "tissue_positions_list.csv"), stringsAsFactors = FALSE,
                          header = FALSE)
  colnames(data_p) <- c("Barcode", "Section", "Spot_Y", "Spot_X",
                      "Image_Y", "Image_X")

  capt_area_p[[i]] <- data_p
}

capt_area <- list()
for (i in sampleNames) {
  data <- read.csv(file.path("./vignettes", sampleDir[i], "outs/spatial",
                                  "tissue_positions_list.csv"), stringsAsFactors = FALSE,
                        header = FALSE)
  colnames(data) <- c("Barcode", "Section", "Spot_Y", "Spot_X",
                                  "Image_Y", "Image_X")

  capt_area[[i]] <- data
}

ggplot() +
  geom_sf(aes(geometry = colGeometry(msfe@sfe_data$JBO019, "spotHex")$geometry)) +
  geom_point(data = capt_area$JBO019,
             aes(x = Image_X, y = Image_Y, colour = as.factor(Section))) +
  labs(title = "Liver JBO019") +
  theme_bw()

ggplot() +
  geom_sf(aes(geometry = colGeometry(msfe_p@sfe_data$H1_4)$geometry)) +
  geom_point(data = capt_area_p$H1_4,
             aes(x = Image_X, y = Image_Y, colour = as.factor(Section))) +
  labs(title = "Prostate H1_4") +
  theme_bw()

#------------------------------------------------------------------------------#
samples <- c(sampleDir[2], sampleDir_p[1])
sample_id <- c(sampleNames[2], sampleNames_p[1])
msfe_t <- list(JBO019 = msfe@sfe_data$JBO019,
               H1_4 = msfe_p@sfe_data$H1_4)


i = 1

res <- "fullres"
res <- STExplorer:::.int_resSwitch(res)

cData <- cbind(colData(msfe_t[[i]]), spatialCoords(msfe_t[[i]]))
# for (i in seq_along(sample_id)) {
  .sp_diam <- S4Vectors::metadata(msfe_t[[i]])$spotDiameter[[sample_id[i]]][[res]]
  data <- read.csv(file.path("./vignettes", samples[i], "outs/spatial",
                             "tissue_positions_list.csv"), stringsAsFactors = FALSE,
                   header = FALSE)
  # rows_with_characters <- apply(data[, -1], 1, function(row) {
  #   any(is.na(sapply(row, .int_char_to_numeric)))
  # })
  # data <- data[!rows_with_characters, ]
  # data[, -1] <- lapply(data[, -1], as.numeric)
  colnames(data) <- c("Barcode", "Section", "Spot_Y", "Spot_X",
                      "Image_Y", "Image_X")
  int_list_minMax <- STExplorer:::.int_spotHex_minMax(data)
  subset_sample <- cData$sample_id == sample_id[i]
  cData_sub <- cData[subset_sample, ]
  int_list_subset <- STExplorer:::.int_spotHex_subset(int_list_minMax,
                                         cData_sub, data)
  int_df_perim <- STExplorer:::.int_spotHex_gen(int_list_subset, .sp_diam)
  dtPerim <- rbind(int_df_perim, data)

  capt_area$JBO019$Image_X_scaled <- 0.01006738*(capt_area$JBO019$Image_X)
  capt_area$JBO019$Image_Y_scaled <- 0.01006738*(capt_area$JBO019$Image_Y)

  ggplot() +
    geom_sf(aes(geometry = colGeometry(msfe@sfe_data$JBO019, "spotHex")$geometry)) +
    geom_point(data = dtPerim,
               aes(x = Image_X, y = Image_Y, colour = as.factor(Section))) +
    labs(title = "Liver JBO019") +
    theme_bw()

  ggplot() +
    geom_point(data = capt_area$JBO019,
               aes(x = Spot_X, y = Spot_Y, colour = as.factor(Section))) +
    labs(title = "Prostate H1_4") +
    theme_bw()

  ggplot() +
    geom_point(data = capt_area$JBO019,
               aes(x = Image_X_scaled, y = Image_Y_scaled, colour = as.factor(Section))) +
    labs(title = "Prostate H1_4") +
    theme_bw()

  i = 2

  res <- "fullres"
  msfe_t[[i]] <- STExplorer:::add.spotCntd(msfe_t[[i]], sample_id = sample_id[i])
  msfe_t[[i]] <- STExplorer:::spot.diameter(sfe = msfe_t[[i]], samples = samples[i], sample_id = sample_id[i], res = res)

  res <- STExplorer:::.int_resSwitch(res)
  cData <- cbind(colData(msfe_t[[i]]), spatialCoords(msfe_t[[i]]))
  # for (i in seq_along(sample_id)) {
  .sp_diam <- S4Vectors::metadata(msfe_t[[i]])$spotDiameter[[sample_id[i]]][[res]]
  print(.sp_diam)
  data <- read.csv(file.path(samples[i], "outs/spatial",
                             "tissue_positions_list.csv"), stringsAsFactors = FALSE,
                   header = FALSE)
  # rows_with_characters <- apply(data[, -1], 1, function(row) {
  #   any(is.na(sapply(row, .int_char_to_numeric)))
  # })
  # data <- data[!rows_with_characters, ]
  # data[, -1] <- lapply(data[, -1], as.numeric)
  colnames(data) <- c("Barcode", "Section", "Spot_Y", "Spot_X",
                      "Image_Y", "Image_X")
  int_list_minMax <- STExplorer:::.int_spotHex_minMax(data)
  subset_sample <- cData$sample_id == sample_id[i]
  cData_sub <- cData[subset_sample, ]
  int_list_subset <- STExplorer:::.int_spotHex_subset(int_list_minMax,
                                                      cData_sub, data)
  int_df_perim <- STExplorer:::.int_spotHex_gen(int_list_subset, .sp_diam)
  dtPerim <- rbind(int_df_perim, data)

  capt_area_p$H1_4$Image_X_scaled <- 0.01006738*(capt_area_p$H1_4$Image_X)
  capt_area_p$H1_4$Image_Y_scaled <- 0.01006738*(capt_area_p$H1_4$Image_Y)

  ggplot() +
    geom_sf(aes(geometry = colGeometry(msfe_p@sfe_data$H1_4)$geometry)) +
    geom_point(data = dtPerim,
               aes(x = Image_X, y = Image_Y, colour = as.factor(Section))) +
    labs(title = "Prostate H1_4") +
    theme_bw()

  ggplot() +
    geom_point(data = capt_area_p$H1_4,
               aes(x = Spot_X, y = Spot_Y, colour = as.factor(Section))) +
    labs(title = "Prostate H1_4") +
    theme_bw()

  ggplot() +
    geom_point(data = capt_area_p$H1_4,
               aes(x = Image_X_scaled, y = Image_Y_scaled, colour = as.factor(Section))) +
    labs(title = "Prostate H1_4") +
    theme_bw()

  # centroids <- as.data.frame(dtPerim) %>% sf::st_as_sf(coords = c("Image_X",
  #                                                             "Image_Y"), remove = TRUE)
  # cntd_union <- sf::st_union(centroids)
  # voronoi <- sf::st_voronoi(cntd_union, bOnlyEdges = TRUE)
  # voronoi_env <- sf::st_intersection(sf::st_cast(voronoi), sf::st_convex_hull(cntd_union))
  # select_cntds <- centroids$Section == 1
  # polygons <- sf::st_polygonize(voronoi_env) %>% sf::st_cast() %>%
  #   sf::st_sf() %>% sf::st_join(., centroids[select_cntds, ],
  #                       join = st_contains, left = FALSE) %>% dplyr::arrange(Barcode) %>%
  #   dplyr::select(geometry)
  # dataList[[i]] <- polygons
# }

