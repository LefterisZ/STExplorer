#' @param gwpca.obj.var is the gwpca object's $var slot: gwpca.obj$var

prop.var <- function(gwpca.obj.var, n.components) {
    return(rowSums(gwpca.obj.var[,1:n.components])/rowSums(gwpca.obj.var))
}
