#' Internal functions for the \code{add.perimeter.slide} function
#' 
#' 

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

.int_spotHex_minMax <- function(.data) {
    int_list <- list(xmax = max(.data$Spot_X),
                     xmin = min(.data$Spot_X),
                     ymax = max(.data$Spot_Y),
                     ymin = min(.data$Spot_Y))
}

## Fetch data and coordinates for the required locations
.int_spotHex_subset <- function(minMaxList, .cData, .data) {
    int_list <- list(bcd_Xmax = 
                         .cData[.cData$array_col == minMaxList[[1]] - 1, "Barcode"],
                     bcd_Xmin = 
                         .cData[.cData$array_col == minMaxList[[2]] - 1, "Barcode"],
                     bcd_Ymax = 
                         .cData[.cData$array_row == minMaxList[[3]] - 1, "Barcode"],
                     bcd_Ymin = 
                         .cData[.cData$array_row == minMaxList[[4]] - 1, "Barcode"])
    
    out_list <- vector("list", length = 4)
    
    for (i in seq_along(int_list)) {
        out_list[[i]] <- .data[.data$Barcode %in% int_list[[i]], ]
    }
    
    names(out_list) <- c("spC_Xmax", "spC_Xmin", "spC_Ymax", "spC_Ymin")
    
    return(out_list)
}

## Generate the perimeter spots only where they are needed
.int_spotHex_gen <- function(subsetList, .sp_diam) {
    which <- !unlist(lapply(subsetList, isEmpty))
    
    subsetListNames <- names(subsetList)
    dtList <- vector("list", length = 4)
    modifiedDataFrames <- list()
    
    for (i in seq_along(subsetListNames)) {
        if (which[subsetListNames[i]]) {
            dtList[[i]] <- as.data.frame(subsetList[[subsetListNames[i]]]) %>% 
                mutate(Barcode = paste0(Barcode, ".perim"), 
                       Section = 0)
            
            if (subsetListNames[i] %in% c("spC_Ymax", "spC_Ymin")) {
                if (subsetListNames[i] == "spC_Ymax") {
                    dtList[[i]]$Image_X <- dtList[[i]]$Image_X - (1.8 * .sp_diam) 
                } else {
                    dtList[[i]]$Image_X <- dtList[[i]]$Image_X + (1.8 * .sp_diam)
                }
            } else {
                if (subsetListNames[i] == "spC_Xmax") {
                    dtList[[i]]$Image_Y <- dtList[[i]]$Image_Y - (1.8 * .sp_diam)
                } else {
                    dtList[[i]]$Image_Y <- dtList[[i]]$Image_Y + (1.8 * .sp_diam)
                }
            }
            
            modifiedDataFrames[[subsetListNames[i]]] <- dtList[[i]]
        }
    }
    
    # Remove NULL elements from the modifiedDataFrames list
    modifiedDataFrames <- modifiedDataFrames[!sapply(modifiedDataFrames, is.null)]
    modifiedDataFrames <- bind_rows(modifiedDataFrames)
    
    return(modifiedDataFrames)
}


