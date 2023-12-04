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
                          sample_id = paste0("sample",
                                             sprintf("%02d",
                                                     seq_along(samples))),
                          res = c("lowres", "hires", "fullres")) {

    ## Set required resolution
    res <- match.arg(res)

    if (!is.null(names(samples))) {
        sample_id <- names(samples)
    }


    for (smpl in seq_along(samples)) {
        ## Import scale factors
        scaleF <- jsonlite::fromJSON(txt = file.path(samples[smpl],
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

        s_id <- sample_id[smpl]

        ## Add info to metadata
        metadata(sfe)$spotDiameter[[s_id]][[name]] <- s_diam
    }

    ## Return
    return(sfe)
}
