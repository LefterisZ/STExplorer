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
        select_cntds <- centroids$Section == 1 | centroids$Section == 0
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
    for (sID in sample_id) {
      selection <- colData(sfe)$sample_id %in% sID
      colGeometry(sfe, sample_id = sID, "spotHex") <- hexes[selection, ]
    }

    return(sfe)

}

