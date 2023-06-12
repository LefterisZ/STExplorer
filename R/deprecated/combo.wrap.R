#' @name combo.wrap
#' @description a GWPCA pipeline helper function to run multile gwpcas and  
#'              and output the results to a list. The list has to be a pre-
#'              defined empty list. This function is essentially a wrapper of 
#'              the \code{gwpca.param.combo()} function. The 
#' @param var.no number of genes to be selected for the PCA.
#' @param k number of PCs to keep
#' @param kernel the kernel to be used as a distance decay weighting function
#' @param list a predifined empty list
#' 
#' @export


combo.wrap <- function(var.no, k, kernel, list) {
    list <- list.append(list, gwpca.param.combo(var.no, k, kernel))
    return(list)
}
