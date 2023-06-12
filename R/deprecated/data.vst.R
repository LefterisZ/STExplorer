#' @name data.vst
#' 
#' @description a wraper around the varianceStabilizingTransformation() function 
#'              from DESeq2. For more information look at the 
#'              varianceStabilizingTransformation {DESeq2} documentation.
#'              
#' @param data A matrix of counts.
#' 
#' @param blind logical, whether to blind the transformation to the experimental
#'              design. blind=TRUE should be used for comparing samples in a 
#'              manner unbiased by prior information on samples, for example to 
#'              perform sample QA (quality assurance). blind=FALSE should be 
#'              used for transforming data for downstream analysis, where the 
#'              full use of the design information should be made. blind=FALSE 
#'              will skip re-estimation of the dispersion trend, if this has 
#'              already been calculated. If many of genes have large differences 
#'              in counts due to the experimental design, it is important to set 
#'              blind=FALSE for downstream analysis.
#' @param fitType in case dispersions have not yet been estimated for object, 
#'                this parameter is passed on to estimateDispersions (options 
#'                described there).
#' 
#' @return a data of vst counts with genes in columns and locations in rows
#' 
#' @export
#' 

data.vst <- function(data, blind = FALSE, fitType = "local"){
    
    # start counting time
    message("Performing a Variance Stabilising Transformation on the count data.
    This might take a while... 
    Grab your self a coffee!")
    start.t <- Sys.time()
    
    # Variance Stabilising Transformation
    vst <- varianceStabilizingTransformation(data, 
                                           blind = blind, 
                                           fitType = fitType)
  
    vst <- as.data.frame(t(vst)) # transpose and transform to df
    
    # stop counting time
    end.t <- Sys.time()
    elapsed.t <- round(difftime(end.t, start.t, units = "mins"), 2)
    message("VST is complete! 
    Elapsed time is: ", elapsed.t, " minutes.")
    
    # return vst
    return(vst)
}
