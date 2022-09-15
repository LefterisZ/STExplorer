gwpca.cv.contrib <- function (x, loc, bw, k = 2, robust = FALSE, kernel = "bisquare", 
                              adaptive = FALSE, p = 2, theta = 0, longlat = F, pcafun = wpca, dMat) 
{
    if (missing(dMat)) 
        DM.given <- F
    else DM.given <- T
    n <- nrow(loc)
    m <- ncol(x)
    w <- array(0, c(n, m, k))
    score <- numeric(n)
    if (robust == FALSE) 
        pcafun = wpca
    else pcafun = rwpca
    for (i in 1:n) {
        if (DM.given) 
            dist.vi <- dMat[, i]
        else {
            dist.vi <- gw.dist(loc, focus = i, p = p, theta = theta, 
                               longlat = longlat)
        }
        wt <- gw.weight(dist.vi, bw, kernel, adaptive)
        wt[i] <- 0
        use <- wt > 0
        wt <- wt[use]
        if (length(wt) <= 1) {
            score[i] <- Inf
            expr <- paste("Too small bandwidth: ", bw)
            warning(paste(expr, "and the CV value can't be given there.", 
                          sep = ", "))
            break
        }
        v <- pcafun(x[use, ], wt, nu = 0, nv = k)$v
        v <- v %*% t(v)
        score[i] <- sum((x[i, ] - x[i, ] %*% v))^2
    }
    score
}
