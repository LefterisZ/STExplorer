#' @description This function calculates the diameter, in pixels, of each Visium 
#' slide spot using the scale factors. The scale factors can be found in the 
#' \code{scalefactors_json.json} output from spaceranger.
#' 
#' It is good practice to input the pathways to the folder where the 
#' \code{scalefactors_json.json} is placed as an object generated using the 
#' \{base} function \code{file.path()}.
#' 
#' @param samples 	a character vector specifying one or more directories, each 
#' corresponding to a 10x Genomics Visium sample (see Details); if provided, 
#' names will be used as sample identifiers.
#' @param sample_id character string specifying unique sample identifiers, one 
#' for each directory specified via samples; ignored if 
#' !is.null(names(samples)).
#' @param obj the SpatialFeaturesExperiment object
#' @param res is the resolution you used to calculate the pixel XY coordinates 
#' for each spot. It can take as values \code{"lowres"} or \code{"hires"} or
#' \code{"fullres"}.
#' 
#' @export

spot_diameter <- function(obj, 
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
        metadata(obj)$spotDiameter[[s_id]][[name]] <- s_diam
    }
    
    ## Return
    return(obj)
}
