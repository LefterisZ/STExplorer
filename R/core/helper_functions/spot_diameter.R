#' @description This function calculates the diameter, in pixels, of each Visium 
#' slide spot using the scale factors. The scale factors can be found in the 
#' \code{scalefactors_json.json} output from spaceranger.
#' 
#' It is good practice to input the pathways to the folder where the 
#' \code{scalefactors_json.json} is placed as an object generated using the 
#' \{base} function \code{file.path()}.
#' 
#' @param sDir the path to the folder where the \code{scalefactors_json.json} is 
#' placed. The value should be an object generated using the \{base} function 
#' \code{file.path()}.
#' @param scale_file a \code{string} with the scale factors file name. Defaults
#' to the default scale factors name from SpaceRnager Output
#' @param res is the resolution you used to calculate the pixel XY coordinates 
#' for each spot. It can take as values either \code{"lowres"} or \code{"hires"}
#' but defaults to \code{"lowres"}.
#' 
#' @export
 

spot_diameter <- function(sDir, scale_file = "scalefactors_json.json", 
                          res = "lowres") {
    ## import scale factors
    scale_f <- jsonlite::fromJSON(txt = file.path(sDir, 
                                                  scale_file))
    
    ## calculate spot diameter
    if (res == "lowres") {
        s_diam <- scale_f$tissue_lowres_scalef * scale_f$spot_diameter_fullres
    } else if (res == "hires") {
        s_diam <- scale_f$tissue_hires_scalef * scale_f$spot_diameter_fullres
    }
    
    ## Return
    return(s_diam)
}
