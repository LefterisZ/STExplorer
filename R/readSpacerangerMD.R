#' @name readSpacerangerMD (read Spaceranger MetaData)
#' @description This function takes the \code{tissue_positions_list.csv} file 
#' from spaceranger and imports it as a data.frame for subsequent use. This 
#' function also adds coordinates as pixels for each slide spot. To do this, the 
#' function, will import the scale factors from the 
#' \code{scale_f_json.json} file which is also a spaceranger output and is
#' located in the same directory as the \code{tissu_positions_list.csv} file.
#' 
#' @param dir The folder path where the tissue positions .csv file is 
#' located.
#' @param file Takes a string. The tissue positions .csv filename. Because the 
#' default spaceranger filename is \code{"tissue_positions_list.csv"} this is 
#' aslo the default value for this parameter.
#' @param res the image resolution you want to continue the analysis with. This 
#' argument is used to calculate the pixel coordinates based on scale factors 
#' for the high resolution or the low resolution image. Takes values either 
#' \code{"high"} or \code{"low"}. Defaults to \code{"low"}. 
#' 
#' 
#' @export

readSpacerangerMD <- function(dir, file = "tissue_positions_list.csv", 
                            res = "low") {
    
    ## read-in the csv file and add colnames
    input_data <- read.csv(file.path(dir, file), 
                           stringsAsFactors = FALSE, header = FALSE)
    colnames(input_data) <- c("Barcode", "Section", "Spot_Y", 
                            "Spot_X", "Image_Y", "Image_X")
    
    ## import the scale factors
    scale_f <- jsonlite::fromJSON(txt = file.path(dir, 
                                                  "scalefactors_json.json"))
    
    ## calculate spot X/Y position in pixels
    if(res == "high"){
        input_data$pixel_x <- input_data$Image_X * scale_f$tissue_hires_scalef
        input_data$pixel_y <- input_data$Image_Y * scale_f$tissue_hires_scalef
    }else if(res == "low"){
        input_data$pixel_x <- input_data$Image_X * scale_f$tissue_lowres_scalef
        input_data$pixel_y <- input_data$Image_Y * scale_f$tissue_lowres_scalef
    }
    
    ## Return
    return(input_data)
}

