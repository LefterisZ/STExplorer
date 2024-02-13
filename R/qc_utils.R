#' Add barcodes in the colData
#'
#' @name add.barcodes
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to add the barcodes in the \code{colData} slot of a
#' \code{SpatialFeaturesExperiment} (SFE) object. In case of multiple samples
#' the \code{read10xVisiumSFE} function adds a suffix in the barcodes to make
#' them unique so that they can be used as col/rownames. This function takes
#' care of this by removing the suffix when transferring the barcodes to a
#' column.
#'
#' @keywords internal
#'
#' @param sfe The SFE object.
#'
#' @details
#' This function adds a new column named "Barcode" to the \code{colData} slot
#' of the given \code{SpatialFeaturesExperiment} (SFE) object. The function
#' moves the existing row names to this new column. If the barcodes in the
#' \code{colData} have a suffix in the format "-(0-9)-(0-9)", which is added
#' by the \code{read10xVisiumSFE} function to make them unique for multiple
#' samples, this function removes the suffix from the barcodes in the new
#' column. Additionally, the function extracts the capture area index from the
#' barcodes and adds it as a new column named "Capt_area" in the \code{colData}.
#'
#' @param sfe The SpatialFeaturesExperiment object.
#'
#' @returns The modified SpatialFeaturesExperiment object with added barcode
#' columns.
#' @importFrom SpatialFeatureExperiment colData
#' @importFrom SpatialFeatureExperiment colData<-
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso
#' \code{\link{read10xVisiumSFE}}
#'
add.barcodes <- function(sfe) {
  ## Move rownames to a column
  colData(sfe)$Barcode <- rownames(colData(sfe))
  ## Find which rows have a suffix
  bc_check <- grepl("-[0-9]-[0-9]$", colData(sfe)$Barcode)
  ## Remove the suffix from these rows in the new column
  colData(sfe)$Barcode[bc_check] <- gsub("-[0-9]$", "",
                                         colData(sfe)$Barcode[bc_check])
  ## Add the capture area index
  colData(sfe)$Capt_area <- gsub("^[A,T,G,C]*-", "",
                                 colData(sfe)$Barcode)

  return(sfe)
}


#' Add ground truth annotation to SpatialFeaturesExperiment object
#'
#' @name add.gTruth
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to add the ground truth annotation in the colData slot of a
#' SpatialFeaturesExperiment (SFE) object. Multiple annotations can be imported
#' using this approach, but note that if not all spots are annotated, NAs will
#' be imported. The SFE object has other slots for adding extra annotations that
#' may not cover the entire dataset.
#'
#' @param sfe The SpatialFeaturesExperiment object.
#' @param gtruth A dataframe containing the ground truth for the dataset.
#' It should have at least three columns: "Barcode" (spot barcodes matching
#' the colnames of the SFE object), "sample_id" (sample ID for each spot,
#' important for multiple samples), and the ground truth annotation column.
#'
#' @return The modified SpatialFeaturesExperiment object with the ground truth
#' annotation added to the colData slot.
#'
#' @details
#' This function adds the ground truth annotation to the colData slot of a
#' SpatialFeaturesExperiment object. It performs a merge between the colData of
#' the SFE object and the provided ground truth dataframe based on "Barcode" and
#' "sample_id". The merge includes all spots, and any missing annotations are
#' represented as NAs in the merged colData. The row names of the colData are
#' temporarily stored in a separate column before merging and then set back as
#' row names after the merge operation.
#'
#' @keywords internal
#'
#' @importFrom SpatialFeatureExperiment colData
#' @importFrom SpatialFeatureExperiment colData<-
#' @importFrom S4Vectors DataFrame
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso
#' \code{\link{add.index}}: A function to add an index to the colData slot of a
#' SpatialFeaturesExperiment object.
#' \code{\link{add.barcodes}}: A function to add barcodes to the colData slot of
#' a SpatialFeaturesExperiment object.
#'
#'
add.gTruth <- function(sfe, gtruth) {
  ## Add the row names in a column because after merging are being lost.
  # colData(sfe)$rowNames <- colnames(sfe) #--> deprecated. Delete?
  ## Merge colData and ground truth
  merger <- merge(colData(sfe), DataFrame(gtruth),
                  by = c("Barcode", "sample_id"),
                  all = TRUE)

  ## Re-order the merger
  ## - The merger is ordered automatically on the combo "Barcode"/"sample_id".
  ##   Thus, when multiple samples are present it creates a missmatch of the
  ##   annotations.
  ## - Here we transform the index column from character to numeric and use it
  ##   to sort the merger as it should have been. If we don't do it then the
  ##   annotations are not matched correctly.
  ## Get the number of unique samples
  n_samples <- length(unique(colData(sfe)[, "sample_id"]))
  if (n_samples > 1) {
    merger$index2 = as.numeric(sub("spot_", "", merger$index))
    merger <- merger[order(merger$index2),]
    merger$index2 <- NULL
  }

  ## Add into colData
  colData(sfe)$annotation <- merger$annotation

  ## Set rownames back
  # rownames(colData(sfe)) <- colData(sfe)$rowNames #--> deprecated. Delete?
  # colData(sfe)$rowNames <- NULL #--> deprecated. Delete?

  return(sfe)
}


#' Add index to SpatialFeaturesExperiment object
#'
#' @name add.index
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to add the index in the colData slot of a
#' SpatialFeaturesExperiment (SFE) object. The index serves as a unique
#' identifier for every spot, even if there are multiple samples in the SFE
#' object.
#'
#' @keywords internal
#'
#' @param sfe The SpatialFeaturesExperiment object.
#'
#' @return The modified SpatialFeaturesExperiment object with the index added to
#' the colData slot.
#'
#' @details
#' This function adds an index column to the colData slot of a
#' SpatialFeaturesExperiment object. The index is generated as "spot_" followed
#' by a numeric sequence, starting from 1 and incrementing for each spot.
#' The index provides a unique identifier for every spot, ensuring each spot can
#' be uniquely identified even in the presence of multiple samples.
#'
#' @importFrom SpatialFeatureExperiment colData
#'
#' @seealso
#' \code{\link{add.barcodes}}: A function to add barcodes to the colData slot of
#' a SpatialFeaturesExperiment object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
add.index <- function(sfe) {
  ## Get the length of the SFE object
  len <- dim(colData(sfe))[1]
  ## Add the index
  colData(sfe)$index <- sprintf("spot_%d", seq(1:len))

  return(sfe)
}


#' Generate centroids from spot coordinates in a SpatialFeaturesExperiment
#' object
#'
#' @name add.spotCntd
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to generate centroids from the spot coordinates. It takes the
#' coordinates out of the \code{spatialCoords} slot of a
#' SpatialFeaturesExperiment (SFE) object and stores it inside the colGeometries
#' slot.
#'
#' @keywords internal
#'
#' @param sfe The SpatialFeaturesExperiment object.
#' @param sample_id A character string specifying unique sample identifiers, one
#' or each directory specified via \code{samples} when you loaded the SFE object
#' using the \code{read10xVisiumSFE} function.
#'
#' @details
#' This function extracts the spot coordinates from the \code{spatialCoords}
#' slot of the SFE object, generates centroids from the coordinates, and stores
#' the resulting centroids in the \code{colGeometries} slot of the SFE object.
#' The centroid coordinates are represented as an \code{sf} object, with each
#' centroid corresponding to a spot in the SFE object. The row names of the
#' \code{cntds} object are set to match the column names of the SFE object for
#' easy reference and linking between the two.
#'
#' @importFrom SpatialFeatureExperiment spatialCoords colGeometry
#' @importFrom SpatialFeatureExperiment colGeometry<-
#' @importFrom sf st_as_sf
#' @importFrom magrittr %>%
#'
#' @seealso
#' \code{\link{add.index}}: A function to add an index in the colData slot of a
#' SpatialFeaturesExperiment object.
#' \code{\link{add.barcodes}}: A function to add barcodes in the colData slot of
#' a SpatialFeaturesExperiment object.
#' \code{\link{add.gTruth}}: A function to add ground truth annotation in the
#' colData slot of a SpatialFeaturesExperiment object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @returns The modified SpatialFeaturesExperiment object with the centroid
#' coordinates stored in the colGeometries slot.
#'
add.spotCntd <- function(sfe, sample_id) {
  ## Fetch coordinates
  coords <- colnames(spatialCoords(sfe))
  ## Get centroids
  cntds <- spatialCoords(sfe) %>%
    as.data.frame() %>%
    st_as_sf(coords = coords,
             remove = TRUE)
  ## Add inside colGeometries slot
  for (sID in sample_id) {
    selection <- colData(sfe)$sample_id %in% sID
    cntd_subset <- cntds[selection, ]
    ## Add rownames
    rownames(cntd_subset) <- colnames(sfe)
    colGeometry(sfe, sample_id = sID, "spotCntd") <- cntd_subset
  }

  return(sfe)
}


#' Add Hexagonal geometries
#'
#' @name add.spotHex
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' This function adds hexagons in the
#' \code{SpatialFeatureExperiment} (SFE) object utilizing a perimeter of
#' off-tissue spots that surround the Visium slide spots. This perimeter assists
#' the Voronoi tessellation process by removing any instances where on-tissue
#' spots at the edges end up without a polygon due to the way the tessellation
#' is calculated. The function also provides options to adjust the perimeter
#' addition based on the tilt of pixel coordinates.
#'
#' @keywords internal
#'
#' @param sfe The \code{SpatialFeaturesExperiment} object.
#' @param samples The samples' folder path as provided at
#' \code{\link{read10xVisiumSFE}}.
#' @param sample_id A character string specifying unique sample identifiers, one
#' or each directory specified via \code{samples} when you loaded the SFE object
#' using the \code{read10xVisiumSFE} function.
#' @param res The desired resolution. Can take one of the following values:
#' "lowres", "hires", "fullres".
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @importFrom dplyr bind_rows arrange select
#' @importFrom sf st_as_sf st_cast st_join st_union st_voronoi st_intersection
#' @importFrom sf st_convex_hull st_polygonize st_sf st_contains
#' @importFrom SpatialFeatureExperiment colGeometry
#' @importFrom S4Vectors metadata
#' @importFrom magrittr %>%
#' @importFrom utils read.csv
#'
#' @seealso \code{\link{read10xVisiumSFE}}
#'

# The idea is to work on the array rows and columns. If there are spots on the
# first and last row/ column then we select those, make a TRUE/FALSE vector to
# fetch the coordinates from the spatialCoords. Then add the diameter*1.5,
# then take their new coordinates and add them back to spatialCoords
# (problems with colData and SFE? Maybe export spatialCoords work and then
# add to SFE ONLY the hexagons.)
add.spotHex <- function(sfe,
                        samples,
                        sample_id,
                        res = c("lowres", "hires", "fullres")) {
  ## Prepare required data
  res <- match.arg(res)
  res <- .int_resSwitch(res)
  cData <- cbind(colData(sfe), spatialCoords(sfe))
  n <- length(sample_id)
  dataList <- vector("list", length = n)

  for (i in seq_along(sample_id)) {
    ## Get spot diameter
    .sp_diam <- metadata(sfe)$spotDiameter[[sample_id[i]]][[res]]
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
    subset_sample <- cData$sample_id == sample_id[i]
    cData_sub <- cData[subset_sample,]
    int_list_subset <- .int_spotHex_subset(int_list_minMax, cData_sub, data)

    ## Generate the perimeter spots
    int_df_perim <- .int_spotHex_gen(int_list_subset, .sp_diam)

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
    select_cntds <- centroids$Section == 1 # centroids$Section == 0
    polygons <- st_polygonize(voronoi_env) %>% # polygonise the tessellation
      st_cast() %>% # convert GEOMETRYCOLLECTION to multiple POLYGONS
      st_sf() %>%  # convert sfc object to sf for st_join afterwards
      st_join(.,
              centroids[select_cntds, ],
              join = st_contains,
              left = FALSE) %>% # Join the centroids with the POLYGONS
      arrange(.data$Barcode) %>%
      dplyr::select(.data$geometry)

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

  return(sfe)

}


#' Calculate gene Coefficient of Variation
#'
#' @name get.QC.CoefficientOfVar
#'
#' @description
#' IMPORTANT: This function is not exported to be used alone.
#' A function to calculate the coefficient of variation (CV) for each gene.
#' CV for each gene is calculated as the ratio of the standard deviation to the
#' mean and is expressed as a percentage. A lower CV would therefore mean higher
#' stability.
#'
#' @keywords internal
#'
#' @param sfe A SpatialFeaturesExperiment object.
#' @param assay The name of the assay to use.
#' @param .sample_id The sample ID to calculate the CV for.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @details This function calculates the CV for each gene in the given
#' `SpatialFeaturesExperiment` object and assigns the CV values to the rowData
#' of the object. The CV is calculated separately for the sample (`s_CV`) and
#' population (`p_CV`) data based on the provided `assay` and `.sample_id`.
#'
#' @return
#' The modified SpatialFeaturesExperiment object with added Coefficient of
#' Variance statistic calculated.
#'
#' @importFrom SpatialFeatureExperiment rowData colData
#' @importFrom SummarizedExperiment assay
#'
#' @seealso \code{\link{SpatialFeatureExperiment}}
#'
get.QC.CoefficientOfVar <- function(sfe, assay, .sample_id) {

  ## Fetch the data per sample and perform the calculations
  data <- assay(sfe, assay)[,colData(sfe)$sample_id == .sample_id]

  ## Calculate CV only for locations the gene is present
  sd = rowData(sfe)[[paste0(.sample_id, ".s_SD")]]
  mean = rowData(sfe)[[paste0(.sample_id, ".s_mean")]]
  s_CV = (sd / mean) * 100

  ## Calculate CV only for all locations in the sample
  sd = rowData(sfe)[[paste0(.sample_id, ".p_SD")]]
  mean = rowData(sfe)[[paste0(.sample_id, ".p_mean")]]
  p_CV = (sd / mean) * 100

  rowData(sfe)[paste0(.sample_id, ".s_CV")] <- s_CV
  rowData(sfe)[paste0(.sample_id, ".p_CV")] <- p_CV

  return(sfe)
}

#' Calculate gene expression statistics
#'
#' @name get.QC.ExprStats
#'
#' @description
#' A function to calculate additional expression statistics.
#'
#' @keywords internal
#'
#' @param sfe A SpatialFeaturesExperiment object.
#' @param assay The name of the assay to use.
#' @param .sample_id The sample ID to run the calculations for.
#'
#' @details
#' This function calculates various expression statistics for a given sample in
#' a SpatialFeaturesExperiment object. The statistics include the number of
#' locations a gene is present, the mean and median UMI counts over the
#' locations, the minimum and maximum UMI counts, and the standard deviation of
#' UMI counts.
#'
#' @return
#' The modified SpatialFeaturesExperiment object with added expression
#' statistics.
#'
#' @importFrom Matrix rowSums
#' @importFrom sparseMatrixStats rowMaxs
#' @importFrom SpatialFeatureExperiment rowData colData
#' @importFrom SummarizedExperiment assay
#' @importFrom stats sd median
#'
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso \code{\link{get.QC.CoefficientOfVar}}
#'
get.QC.ExprStats <- function(sfe, assay, .sample_id) {

  ## Fetch the data per sample and perform the calculations
  data <- assay(sfe, assay)[,colData(sfe)$sample_id == .sample_id]

  ## Get the number of locations a gene is present
  s_nLocations <- Matrix::rowSums(data != 0)
  p_nLocations <- sum(colData(sfe)$sample_id == .sample_id)

  ## Get the mean of UMI counts over the locations a gene is present
  s_mean <- rowData(sfe)[[paste0(.sample_id, ".total")]] / s_nLocations
  ## Get the mean of UMI counts over the whole sample
  p_mean <- rowData(sfe)[[paste0(.sample_id, ".total")]] / p_nLocations

  ## Get the max of UMI counts of the locations a gene is present
  max <- sparseMatrixStats::rowMaxs(data)

  ## Get the median of UMI counts over the locations a gene is present
  data <- Matrix::as.matrix(data)
  s_median <- apply(data, 1,
                    function(dt) {
                      dt <- dt[dt != 0]
                      median <- median(dt)
                    })
  ## Get the median of UMI counts over the whole sample
  p_median <- apply(data, 1, median)

  ## Get the min of UMI counts of the locations a gene is present
  s_min <- apply(data, 1,
                 function(dt) {
                   dt <- dt[dt != 0]
                   min <- min(dt)
                 })

  ## Get the Stand. Dev. of UMI counts of the locations a gene is present
  s_SD <- apply(data, 1,
                function(dt) {
                  dt <- dt[dt != 0]
                  sd <- sd(dt)
                })

  ## Get the Stand. Dev. of UMI counts of the whole sample
  p_SD <- apply(data, 1, sd)

  ## Add them to the rowData
  rowData(sfe)[paste0(.sample_id, ".nLocations")] <- s_nLocations
  rowData(sfe)[paste0(.sample_id, ".s_min")] <- s_min
  rowData(sfe)[paste0(.sample_id, ".max")] <- max
  rowData(sfe)[paste0(.sample_id, ".s_mean")] <- s_mean
  rowData(sfe)[paste0(.sample_id, ".s_median")] <- s_median
  rowData(sfe)[paste0(.sample_id, ".s_SD")] <- s_SD
  rowData(sfe)[paste0(.sample_id, ".p_mean")] <- p_mean
  rowData(sfe)[paste0(.sample_id, ".p_median")] <- p_median
  rowData(sfe)[paste0(.sample_id, ".p_SD")] <- p_SD

  return(sfe)
}


#' Find non-expressed genes
#'
#' @name get.QC.FindZeroExpr
#'
#' @description
#' A function to find zero expression values and calculate total UMIs.
#'
#' @keywords internal
#'
#' @param sfe A SpatialFeatureExperiment object.
#'
#' @param assay The name of the assay to use.
#'
#' @param .sample_id The sample ID to run the calculations for.
#'
#' @details This function takes a SpatialFeaturesExperiment object and
#' calculates the total number of UMIs for each gene in the specified sample.
#' It extracts the data for the given sample from the specified assay, finds
#' the total UMIs for each gene, and stores the results in the rowData of the
#' SpatialFeaturesExperiment object.
#'
#' @importFrom Matrix rowSums
#' @importFrom SpatialFeatureExperiment rowData
#' @importFrom SummarizedExperiment assay
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso \code{\link{SpatialFeatureExperiment}}, \code{\link{rowData}},
#' \code{\link{assay}}
#'
#' @return The modified SpatialFeaturesExperiment object with the total UMIs
#' stored in the rowData.
#'

get.QC.FindZeroExpr <- function(sfe, assay, .sample_id){

  ## Fetch the data per sample and perform the calculations
  data <- assay(sfe, assay)[,colData(sfe)$sample_id == .sample_id]

  ## Find the total in-sample number of UMIs
  total <- Matrix::rowSums(data)

  rowData(sfe)[paste0(.sample_id, ".total")] <- total

  return(sfe)
}


#' Calculate sparsity
#'
#' @name get.QC.Sparsity
#'
#' @description
#' A function to calculate the sparsity of the dataset and get the
#' sparsity of each location or each gene.
#'
#' @keywords internal
#'
#' @param sfe A SpatialFeaturesExperiment object.
#'
#' @param assay The name of the assay to use.
#'
#' @param MARGIN Specifies the aspect for which the sparsity will be calculated.
#' Use 1 for features (genes) or 2 for locations.
#'
#' @param sampleNo The sample number used to call the sparsity over all datasets
#' when multiple samples exist. This parameter is used only when the function
#' is called by addPerGeneQC.
#'
#' @param .sample_id The sample ID to run the calculations for. This parameter
#' is used only when the function is called by addPerGeneQC.
#'
#' @details The sparsity is calculated as the proportion of zeros in the dataset.
#' It provides an indication of the density of expression values. A higher
#' sparsity indicates a higher number of zero values, indicating low expression
#' levels or absence of expression.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @seealso addPerGeneQC
#'
get.QC.Sparsity <- function(sfe,
                            assay,
                            MARGIN,
                            sampleNo = NULL,
                            .sample_id = NULL) {

  ## Fetch data per sample
  if (MARGIN == 1) {
    sfe <- .int_sparsity_1(sfe = sfe, assay = assay,
                           sampleNo = sampleNo, .sample_id = .sample_id)
  } else if (MARGIN == 2) {
    sfe <- .int_sparsity_2(sfe = sfe, assay = assay)
  }

  return(sfe)
}


#' Calculate spot diameter
#'
#' @name spot.diameter
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' It is used internally to calculate the diameter, in pixels, of each Visium
#' slide spot using the scale factors. The scale factors can be found in the
#' "scalefactors_json.json" output from spaceranger. It is recommended to
#' provide the pathways to the folder where the "scalefactors_json.json" files
#' are located using the samples parameter generated using the file.path
#' function.
#'
#' @param sfe The SpatialFeaturesExperiment object.
#' @param samples A character vector specifying one or more directories, each
#' corresponding to a 10x Genomics Visium sample. If provided, the names will
#' be used as sample identifiers.
#' @param sample_id A character string specifying unique sample identifiers,
#' one for each directory specified via samples. This parameter is ignored if
#' names(samples) is not NULL.
#' @param res The resolution used to calculate the pixel XY coordinates for
#' each spot. It can take the values "lowres", "hires", or "fullres".
#'
#' @keywords internal
#'
#' @importFrom jsonlite fromJSON
#'
#' @details
#' This function calculates the diameter, in pixels, of each Visium slide spot
#' using the scale factors provided in the "scalefactors_json.json" files.
#' The scale factors are imported from the specified directories (samples). The
#' res parameter specifies the resolution to use for the calculations:
#' "lowres", "hires", or "fullres". The calculated spot diameters are added to
#' the metadata of the sfe object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
spot.diameter <- function(sfe,
                          samples,
                          sample_id,
                          res = c("lowres", "hires", "fullres")) {

  ## Set required resolution
  res <- match.arg(res)

  # if (!is.null(names(samples))) {
  #   sample_id <- names(samples)
  # }

  #for (smpl in seq_along(samples)) {
  ## Import scale factors
  scaleF <- jsonlite::fromJSON(txt = file.path(samples[sample_id],
                                               "outs/spatial",
                                               "scalefactors_json.json"))

  ## Calculate spot diameter
  if (res == "lowres") {
    s_diam <- scaleF$tissue_lowres_scalef * scaleF$spot_diameter_fullres
    name <- "spot_diameter_lowres"
  } else if (res == "hires") {
    s_diam <- scaleF$tissue_hires_scalef * scaleF$spot_diameter_fullres
    name <- "spot_diameter_hires"
  } else if (res == "fullres") {
    s_diam <- scaleF$spot_diameter_fullres
    name <- "spot_diameter_fullres"
  }

  #s_id <- sample_id[smpl]

  ## Add info to metadata
  metadata(sfe)$spotDiameter[[sample_id]][[name]] <- s_diam
  #}

  ## Return
  return(sfe)
}


# ---------------------------------------------------------------------------- #
#  ################ INTERNAL FUNCTIONS ASSOCIATED RESOLUTION ################
# ---------------------------------------------------------------------------- #
#' Internal: switch resolution name
#'
#' @name dot-int_resSwitch
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function that selects the name of the spot diameter from the required
#' resolution which is stored inside SpatialFeatureExperiment's metadata.
#'
#' @keywords internal
#'
#' @param .res resolution name
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
.int_resSwitch <- function(.res) {
  if (.res == "lowres") {
    name <- "spot_diameter_lowres"
  } else if (.res == "hires") {
    name <- "spot_diameter_hires"
  } else if (.res == "fullres") {
    name <- "spot_diameter_fullres"
  }

  return(name)
}


# ---------------------------------------------------------------------------- #
#  ########### INTERNAL FUNCTIONS ASSOCIATED WITH HEXAGON GEOMS #############
# ---------------------------------------------------------------------------- #
## Internal functions for the \code{add.perimeter.slide} function
##
##

## ## Get min/max values for the first/last capture area rows/columns. It
## is not fixed because Visium can have 6.5mm or 11mm capture areas. The 6.5mm
## has column indexes that range from 0 to 127 and row indexes that range from
## 0 to 77. The 11mm has column indexes that range from 0 to 223 and row indexes
## that range from 0 to 128.

## ## It has come to our attention that when loading an experiment with the
## `SpatialFeatureExperiment`'s `read10XVisiumSFE` function, the numbering of
## the array cols and rows starts from 1 instead of 0. As a result, the maximums
## that we mention above are actually maximum - 1.

## ## Another point is that the image is rotated 90 degrees and flipped. Such as
## the array columns are X coordinates in pixels and the array rows are the Y
## coordinates in pixels and the min and max are inverted.

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Internal: get min and max array col/row
#'
#' @name dot-int_spotHex_minMax
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' This function calculates the minimum and maximum values of x and y
#' coordinates for the spots in the provided data.
#'
#' @param .data A data frame containing the spot coordinates.
#'
#' @returns A list with the following elements:
#' \describe{
#' \item{xmax}{The maximum x-coordinate value.}
#' \item{xmin}{The minimum x-coordinate value.}
#' \item{ymax}{The maximum y-coordinate value.}
#' \item{ymin}{The minimum y-coordinate value.}
#' }
#'
#' @details
#' This function is used internally to calculate the minimum and maximum x and y
#' coordinates for the spots in the provided data frame.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
.int_spotHex_minMax <- function(.data) {
  int_list <- list(xmax = max(.data$Spot_X),
                   xmin = min(.data$Spot_X),
                   ymax = max(.data$Spot_Y),
                   ymin = min(.data$Spot_Y))
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Internal: get info for spots on the edge
#'
#' @name dot-int_spotHex_subset
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' This function fetches data and coordinates for the required locations
#' based on the provided minimum and maximum values.
#'
#' @param minMaxList A list containing the minimum and maximum x and y
#' coordinate values.
#'
#' @param .cData A data frame containing the barcode and array coordinate
#' information.
#'
#' @param .data A data frame containing the spot data.
#'
#' @importFrom dplyr filter
#'
#' @returns A list with the following elements:
#' \describe{
#' \item{bcd_Xmax}{Data frame with the spots corresponding to the maximum
#' x-coordinate value.}
#' \item{bcd_Xmin}{Data frame with the spots corresponding to the minimum
#' x-coordinate value.}
#' \item{bcd_Ymax}{Data frame with the spots corresponding to the maximum
#' y-coordinate value.}
#' \item{bcd_Ymin}{Data frame with the spots corresponding to the minimum
#' y-coordinate value.}
#' \item{bcd_Xmax1}{Data frame with the spots corresponding to the x-coordinate
#' value adjacent to the maximum x-coordinate.}
#' \item{bcd_Xmin1}{Data frame with the spots corresponding to the x-coordinate
#' value adjacent to the minimum x-coordinate.}
#' }
#'
#' @details
#' This function is used internally to fetch data and coordinates for specific
#' locations based on the provided minimum and maximum values. It takes the
#' barcode and array coordinate information from the \code{.cData} data frame
#' and the spot data from the \code{.data} data frame. It returns a list of data
#' frames, where each data frame corresponds to a specific location based on the
#' minimum and maximum values.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
## Fetch data and coordinates for the required locations
.int_spotHex_subset <- function(minMaxList, .cData, .data) {
  int_list <- list(bcd_Xmax =
                     .cData[.cData$array_col == minMaxList[[1]], "Barcode"],
                   bcd_Xmin =
                     .cData[.cData$array_col == minMaxList[[2]], "Barcode"],
                   bcd_Ymax =
                     .cData[.cData$array_row == minMaxList[[3]], "Barcode"],
                   bcd_Ymin =
                     .cData[.cData$array_row == minMaxList[[4]], "Barcode"],
                   bcd_Xmax1 =
                     .cData[.cData$array_col == minMaxList[[1]] - 1, "Barcode"],
                   bcd_Xmin1 =
                     .cData[.cData$array_col == minMaxList[[2]] + 1, "Barcode"],
                   bcd_Ymax1 =
                     .cData[.cData$array_row == minMaxList[[3]] - 1, "Barcode"],
                   bcd_Ymin1 =
                     .cData[.cData$array_row == minMaxList[[4]] + 1, "Barcode"])

  out_list <- vector("list", length = 8)

  for (i in seq_along(int_list)) {
    out_list[[i]] <- .data[.data$Barcode %in% int_list[[i]], ]
  }

  names(out_list) <- c("spC_Xmax", "spC_Xmin", "spC_Ymax",
                       "spC_Ymin", "spC_Xmax1", "spC_Xmin1",
                       "spC_Ymax1", "spC_Ymin1")

  return(out_list)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Internal: generate perimeter spots
#'
#' @name int_spotHex_gen
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' This function generates the perimeter spots only where they are needed based
#' on the provided subset list and spot diameter.
#'
#' @keywords internal
#'
#' @param subsetList A list containing subsets of spot data frames.
#'
#' @param .sp_diam The diameter of the spots.
#'
#' @importFrom dplyr mutate bind_rows
#' @importFrom magrittr %>%
#'
#' @returns A data frame containing the modified spot data frames for the
#' perimeter spots.
#'
#' @details
#' This function is used internally to generate the perimeter spots only for the
#' locations where they are needed. It takes a subset list containing subsets of
#' spot data frames and modifies them to create the perimeter spots. The
#' modification includes adding a ".perim" suffix to the Barcode column, setting
#' the Section value to 0, and adjusting the Image_X and Image_Y coordinates
#' based on the subset. The modified spot data frames are then combined into a
#' single data frame and returned.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
## Generate the perimeter spots only where they are needed
.int_spotHex_gen <- function(subsetList, .sp_diam) {
  which <- !unlist(lapply(subsetList, isEmpty))

  subsetListNames <- names(subsetList)
  dtList <- vector("list", length = 8)
  modifiedDFs <- list()

  for (i in seq_along(subsetListNames)) {
    name <- subsetListNames[i]
    if (which[subsetListNames[i]]) {
      dtList[[i]] <- as.data.frame(subsetList[[name]]) %>%
        mutate(Barcode = paste0(.data$Barcode, ".perim"),
               Section = 2)

      if (name %in% c("spC_Ymax", "spC_Ymax1", "spC_Ymin", "spC_Ymin1")) {
        if (name %in% c("spC_Ymax", "spC_Ymax1")) {
          dtList[[i]]$Image_Y <- dtList[[i]]$Image_Y + (3 * .sp_diam)
        } else {
          dtList[[i]]$Image_Y <- dtList[[i]]$Image_Y - (3 * .sp_diam)
        }
      } else {
        if (name %in% c("spC_Xmax", "spC_Xmax1")) {
          dtList[[i]]$Image_X <- dtList[[i]]$Image_X + (1.8 * .sp_diam)
        } else {
          dtList[[i]]$Image_X <- dtList[[i]]$Image_X - (1.8 * .sp_diam)
        }
      }

      modifiedDFs[[name]] <- dtList[[i]]
    }
  }

  # Remove NULL elements from the modifiedDataFrames list
  modifiedDFs <- modifiedDFs[!sapply(modifiedDFs, is.null)]
  modifiedDFs <- bind_rows(modifiedDFs)

  return(modifiedDFs)
}


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#' Internal: Convert a Character to a Numeric Value
#'
#' @name dot-int_char_to_numeric
#'
#' @description
#' IMPORTANT: This function takes a character input and attempts to convert it
#' into a numeric value. If the conversion is successful, it returns the numeric
#' value; otherwise, it returns NA.
#'
#' @keywords internal
#'
#' @param x A character value to be converted to numeric.
#'
#' @returns If the conversion is successful, a numeric value is returned. If
#' the input cannot be converted, NA (Not Available) is returned.
#'
#' @examples
#' char_to_numeric("123")  # Returns: 123
#' char_to_numeric("abc")  # Returns: NA
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
## Function to convert a character to a numeric value,
##  returning NA if it's not a number.
.int_char_to_numeric <- function(x) {
  if (is.na(as.numeric(x))) {
    return(NA)
  } else {
    return(as.numeric(x))
  }
}

# ---------------------------------------------------------------------------- #
#  ############## INTERNAL FUNCTIONS ASSOCIATED WITH SPARSITY ###############
# ---------------------------------------------------------------------------- #
#' Fetch global sparsity stats
#'
#' @name .get.Global.Sparsity
#' @description
#' A function to fetch the dataset sparsity attributes from an assay.
#'
#' @keywords internal
#'
#' @param sfe A SpatialFeatureExperiment object
#'
#' @importFrom S4Vectors metadata
#' @importFrom dplyr bind_rows
#'
.get.Global.Sparsity <- function(sfe) {

  tmp <- metadata(sfe)$Sparsity

  tmp <- dplyr::bind_rows(tmp)

  colnames(tmp) <- c("Assay", "Sample ID", "Total", "Zeros", "Sparsity %")

  return(tmp)
}


#' Internal: calculate feature sparsity
#'
#' @name dot-int_sparsity_1
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to calculate the sparsity of the dataset and get the
#' sparsity of each gene.
#'
#' @keywords internal
#'
#' @param sfe a SpatialFeaturesExperiment objects.
#'
#' @param assay the name of the assay to use.
#'
#' @param sampleNo used to call the sparsity also over aall datasets when
#' multiple samples exist.Used only when MARGIN = 1 at \code{addPerGeneQC}
#' function.
#'
#' @param .sample_id the sample id to run the calculations for. Used only when
#' MARGIN = 1 at \code{addPerGeneQC} function.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialFeatureExperiment colData rowData
#' @importFrom S4Vectors metadata metadata<-
#'
.int_sparsity_1 <- function(sfe,
                            assay,
                            sampleNo,
                            .sample_id) {
  ## Fetch data per sample
  data <- assay(sfe, assay)[,colData(sfe)$sample_id == .sample_id]

  ## Calculate sparsity
  ## per feature over all samples
  if (sampleNo > 1) {
    tot_zeros <- rowSums(assay(sfe, assay) == 0)
    tot_sparsity <- tot_zeros/ncol(assay(sfe, assay))
  }
  ## per feature per sample
  zeros <- rowSums(data == 0)
  sparsity <- zeros/ncol(data)

  ## Add the result to the output matrix
  ## per feature
  if (sampleNo > 1) {
    rowData(sfe)$sparsity_tot <- tot_sparsity
  }
  ## per feature per sample
  rowData(sfe)[paste0(.sample_id, ".sparsity")] <- sparsity

  ## Calculate dataset's sparsity of each location (per sample)
  zeros <- sum(data == 0)
  total <- dim(data)[1] * dim(data)[2]
  sparsity <- round(zeros/total*100, 2)
  metaSpar <- data.frame(assay, .sample_id, zeros, total, sparsity)

  ## Add dataset-wide, per locations, sparsity stats
  metadata(sfe)$Sparsity[[paste0(assay, .sample_id)]] <- metaSpar

  return(sfe)
}


#' Internal: calculate location sparsity
#'
#' @name dot-int_sparsity_2
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function to calculate the sparsity of the dataset and get the
#' sparsity of each location.
#'
#' @keywords internal
#'
#' @param sfe a SpatialFeaturesExperiment objects.
#'
#' @param assay the name of the assay to use.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialFeatureExperiment colData
#' @importFrom Matrix colSums
#'
.int_sparsity_2 <- function(sfe, assay) {
  ## Fetch data per sample
  data <- assay(sfe, assay)

  ## Calculate sparsity
  ## per location
  zeros <- Matrix::colSums(data == 0)
  sparsity <- zeros/nrow(data)


  ## Add the result to the output matrix
  ## per location
  colData(sfe)$sparsity <- sparsity

  return(sfe)
}
