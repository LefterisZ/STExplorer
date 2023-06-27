#' Internal: switch resolution name
#'
#' @name .int_resSwitch
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' A function that selects the name of the spot diameter from the required
#' resolution which is stored inside SpatialFeatureExperiment's metadata.
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
