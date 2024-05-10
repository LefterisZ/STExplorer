# Test addGeometries
dir_test <- "./data/test_data/Visium_Human_Prostate/Patient1/V1_2"
names(dir_test) <- "V1_2"
res = "lowres"
samples = dir_test
sample_id = "V1_2"
flipped = FALSE
barcodes = "input"
i = 1


## Add Centroids
sfe <- add.spotCntd(sfe,
                    sample_id = sample_id)

## Get/ calculate spot diameter
sfe <- spot.diameter(sfe = sfe,
                     samples = samples,
                     sample_id = sample_id,
                     res = res)

## Prepare required data
# res <- match.arg(res)
res <- .int_resSwitch(res)
cData <- cbind(colData(sfe), spatialCoords(sfe))
n <- length(sample_id)
dataList <- vector("list", length = n)

for (i in seq_along(sample_id)) {
  ## Get spot diameter
  sp_diam <- metadata(sfe)$spotDiameter[[sample_id[i]]][[res]]
  data <- read.csv(file.path(samples[i], "outs/spatial",
                             "tissue_positions_list.csv"),
                   stringsAsFactors = FALSE, header = FALSE)
  ## Check if for whatever reason there is still column names in the first row
  ## Check for rows with characters (non-numeric values) in all columns but the
  ## first
  rows_with_characters <- apply(data[, -1], 1,
                                function(row){
                                  any(is.na(sapply(row,
                                                   .int_char_to_numeric)))
                                }
  )

  ## Remove rows with characters (non-numeric values) at any column (but the
  ## first). We are left with rows that have all columns containing numeric
  ## values and the first column containing the spot barcodes
  data <- data[!rows_with_characters, ]

  ## Convert all numeric columns (except the first) to numeric
  data[, -1] <- lapply(data[, -1], as.numeric)

  ## Add colnames
  colnames(data) <- c("Barcode", "Section", "Spot_Y",
                      "Spot_X", "Image_Y", "Image_X")
  ## Get min/max values for the first/last capture area rows/columns. It
  ## is not fixed because Visium can have 6.5mm or 11mm capture areas.
  int_list_minMax <- .int_spotHex_minMax(data)

  ## Fetch data and coordinates for the required locations
  # # Check if array_col and array_row are not present in cData
  # if (!("array_col" %in% colnames(cData)) &
  #     !("array_row" %in% colnames(cData))) {
  #   # Add Spot_X and Spot_Y columns from data and rename them
  #   cData$array_col <- data$Spot_X
  #   cData$array_row <- data$Spot_Y
  # }

  subset_sample <- cData$sample_id == sample_id[i]
  cData_sub <- cData[subset_sample,]
  int_list_subset <- .int_spotHex_subset(int_list_minMax, cData_sub, data)

  ## Generate the perimeter spots
  int_df_perim <- .int_spotHex_gen(int_list_subset,
                                   sp_diam,
                                   flipped = flipped)

  ## Add them to the rest of the data
  dtPerim <- rbind(int_df_perim,
                   data)

  ## Convert spots to centroids
  centroids <- as.data.frame(dtPerim) %>%
    st_as_sf(coords = c("Image_X", "Image_Y"),
             remove = TRUE)

  ## Combine the points into a multipoint geometry:
  cntd_union <- st_union(centroids)

  ## Use the union of points to generate a voronoi object
  voronoi <- st_voronoi(cntd_union, bOnlyEdges = TRUE)

  ## Create an enveloped voronoi tessellation around the tissue
  voronoi_env <- st_intersection(st_cast(voronoi),
                                 st_convex_hull(cntd_union))

  ## Generate the POLYGONS from the MULTILINESTRING
  # First select the barcodes for which the polygons will be generated
  if (barcodes == "all") {
    select_cntds <- centroids$Section == 1 # centroids$Section == 0
  } else if (barcodes == "input") {
    select_cntds <- centroids$Barcode %in% cData$Barcode
  } else if (length(barcodes) > 1) {
    select_cntds <- centroids$Barcode %in% barcodes
  }
  polygons <- st_polygonize(voronoi_env) %>% # polygonise the tessellation
    st_cast() %>% # convert GEOMETRYCOLLECTION to multiple POLYGONS
    st_sf() %>%  # convert sfc object to sf for st_join afterwards
    st_join(.,
            centroids[select_cntds, ],
            join = st_contains,
            left = FALSE) %>% # Join the centroids with the POLYGONS
    arrange(Barcode) %>%
    dplyr::select(geometry)

  ## Append it to a list
  dataList[[i]] <- polygons
}

hexes <- dplyr::bind_rows(dataList)
rownames(hexes) <- colnames(sfe)
## Add inside colGeometries slot
## NOTE: remove this for loop. Not needed anymore
for (sID in sample_id) {
  selection <- colData(sfe)$sample_id %in% sID
  colGeometry(sfe, sample_id = sID, "spotHex") <- hexes #hexes[selection, ]
}


################################################################################
# Clear environment ------------------------------------------------------------
rm(res, samples, sample_id, flipped, cData, n, dataList, sp_diam, data,
   rows_with_characters, int_list_minMax, subset_sample, cData_sub, barcodes,
   int_list_subset, int_df_perim, dtPerim, centroids, cntd_union, voronoi,
   voronoi_env, select_cntds, polygons, hexes, selection, i, sID, dir_test,
   cntds, selection, cntd_subset, id, merger)
