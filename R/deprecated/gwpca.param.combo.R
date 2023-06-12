#' @name gwpca.param.combo
#' @description a GWPCA pipeline helper function to run gwpca from a selection 
#'              of variables from a table of variable combinations. This table 
#'              can be generated with the \code{param.combo()} function.
#' @param var.no number of genes to be selected for the PCA.
#' @param k number of PCs to keep
#' @param kernel the kernel to be used as a distance decay weighting function
#' 
#' @export

gwpca.param.combo <- function(var.no, k, kernel){
    inputPCAgw <- SpatialPointsDataFrame(coords, vst_df, match.ID = TRUE)
    print(var.no)
    select <- order(row_vars, decreasing = TRUE)[seq_len(var.no)]
    inputPCAgw <- inputPCAgw[select]
    vars <- colnames(inputPCAgw@data)
    bw <- 6*spot_diameter(spatialDir)
    k <- k
    obj <- gwpca(inputPCAgw, 
                 vars = vars, 
                 bw = bw,
                 k = k,
                 kernel = kernel)
    
    return(obj)
}
