#' @name gwpca.plot.ex.outlier
#' @description A function to visualise the genes that make the spot in question
#'              have a high discrepancy score.
#' @param data  a Spatial*DataFrame
#' @param i     a vector of variable names to be evaluated
#' @param loc   the centroid coordinates to calculate distances
#' @param bw    bandwidth used in the weighting function
#' @param ylim  the y limits of the plot
#' @param ylab  a label for the y axis
#' @param fixtrans if TRUE, the transparency of the neighbouring observation 
#'                 plot lines increases with distance; If FALSE a standard 
#'                 (non-spatial) parallel coordinate plot is returned.
#' @param ...   other graphical parameters, (see par)
#' 
#' @export


gwpca.plot.ex.outlier <- function(data, i, loc, bw, ylim = NULL, ylab = "", 
                                  fixtrans = FALSE, ...) {
    m <- ncol(data)
    #bw <- bw*bw
    bw <- bw
    #dists <- rowSums(sweep(loc,2,loc[i,])**2)
    dists <- rowSums(sweep(loc, 2, loc[i,]))
    #nbrlist <- which(dists < bw*bw)
    nbrlist <- which(dists < bw)
    nbrlist <- nbrlist[nbrlist != i]
    wts <- (1 - dists/(bw*bw))^12
    xss <- scale(data)
    span <- 1:m
    tsc <- 50/length(nbrlist)
    
    
    if (is.null(ylim)) ylim <- c(min(xss[nbrlist,]),max(xss[nbrlist,]))
    
    plot(span,xss[i,],type='l',ylim=ylim,
         xlim=c(0.5,m+0.5),col='red',lwd=6,axes=FALSE,xlab="",ylab=ylab)
    #axis(1,at=1:m,labels=colnames(data),las=2,cex.axis=1.2)
    axis(2,at=seq(floor(ylim[1]),ceiling(ylim[2]),by=1),cex.axis=1.2)
    #abline(v=1:m,col=grey(0.6))
    lines(c(1,m),c(0,0),col=grey(0.6))
    
    if (fixtrans) {
        for (nbr in nbrlist) lines(span,xss[nbr,],col=rgb(0,0,0,0.3), lwd=3) }
    else {
        for (nbr in nbrlist) lines(span,xss[nbr,],col=rgb(0,0,0,tsc*wts[nbr]), lwd=3)} 
    
}


data = data.mat
bw = 3*spot_diameter(spatialDir)
i = which(discrepancy_df == max(discrepancy_df))
st_geometry(polygons) <- "geom_cntd"
loc = st_coordinates(polygons)

