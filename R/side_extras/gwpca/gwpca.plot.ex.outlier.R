gwpca.plot.ex.outlier <- function(x, i, loc, bw, ylim = NULL, ylab = "", 
                                  fixtrans = FALSE, ...) {
    m <- ncol(x)
    bw <- bw*bw
    dists <- rowSums(sweep(loc,2,loc[i,])**2)
    nbrlist <- which(dists < bw*bw)
    wts <- (1 - dists/(bw*bw))^12
    xss <- scale(x)
    span <- 1:m
    tsc <- 25/length(nbrlist)
    
    
    if (is.null(ylim)) ylim <- c(min(xss[nbrlist,]),max(xss[nbrlist,]))
    
    plot(span,xss[i,],type='l',ylim=ylim,
         xlim=c(0.5,m+0.5),col='red',lwd=6,axes=FALSE,xlab="",ylab=ylab,...)
    axis(1,at=1:m,labels=colnames(x),las=2,cex.axis=1.2)
    axis(2,at=seq(floor(ylim[1]),ceiling(ylim[2]),by=1),cex.axis=1.2)
    abline(v=1:m,col=grey(0.6))
    lines(c(1,m),c(0,0),col=grey(0.6))
    
    if (fixtrans) {
        for (nbr in nbrlist) lines(span,xss[nbr,],col=rgb(0,0,0,0.3), lwd=3) }
    else {
        for (nbr in nbrlist) lines(span,xss[nbr,],col=rgb(0,0,0,tsc*wts[nbr]), lwd=3)} 
    
}


x = data.mat
