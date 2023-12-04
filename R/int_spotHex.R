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
