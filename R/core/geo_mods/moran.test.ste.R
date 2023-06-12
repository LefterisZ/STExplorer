#' @moran.test.ste
#' 
#' @description A modification of the moran.test function from spdep
#' 
#' @export

moran.test.ste <- function(x, 
                       listw, 
                       randomisation = TRUE, 
                       zero.policy = NULL,
                       alternative = "greater", 
                       rank  =  FALSE, 
                       na.action = na.fail, 
                       spChk = NULL, 
                       adjust.n = TRUE, 
                       drop.EI2 = FALSE) {
    
    alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
    
    if (!inherits(listw, "listw")) {
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    } 
    if (!is.numeric(x)) {
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    }
    
    if (is.null(zero.policy)) {
        zero.policy <- get("zeroPolicy", envir  =  spdep:::.spdepOptions)
    }
    
    stopifnot(is.logical(zero.policy))
    
    if (is.null(spChk)) {
        spChk <- spdep::get.spChkOption()
    }
    
    if (spChk && !chkIDs(x, listw)) {
        stop("Check of data and weights ID integrity failed")
    }
    
    xname <- deparse(substitute(x))
    wname <- deparse(substitute(listw))
    NAOK <- deparse(substitute(na.action)) == "na.pass"
    x <- na.action(x)
    na.act <- attr(x, "na.action")
    
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    
    n <- length(listw$neighbours)
    
    if (n != length(x)) {
        stop("objects of different length")
    } 
    
    wc <- spdep::spweights.constants(listw, zero.policy = zero.policy, 
                              adjust.n = adjust.n)
    S02 <- wc$S0*wc$S0
    res <- spdep::moran(x, listw, wc$n, wc$S0, zero.policy = zero.policy, 
                 NAOK = NAOK)
    I <- res$I
    K <- res$K
    
    if (rank) {
        K <- (3*(3*wc$n^2 - 7))/(5*(wc$n^2 - 1))
    } 
    
    EI <- (-1) / wc$n1
    
    if (randomisation) {
        VI <- wc$n*(wc$S1*(wc$nn - 3*wc$n + 3) - wc$n*wc$S2 + 3*S02)
        tmp <- K*(wc$S1*(wc$nn - wc$n) - 2*wc$n*wc$S2 + 6*S02)
        
        if (tmp > VI) {
            warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
        }
        
        VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
        
        if (!drop.EI2) {
            VI <- (VI - EI^2)
        } 
        
        if (VI < 0) {
            warning("Negative variance,\ndistribution of variable does not meet test assumptions")
        } 
    
    } else {
        VI <- (wc$nn*wc$S1 - wc$n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
        
        if (!drop.EI2) {
            VI <- (VI - EI^2)
        } 
        if (VI < 0) {
            warning("Negative variance,\ndistribution of variable does not meet test assumptions")
        }
    }
    
    ZI <- (I - EI) / sqrt(VI)
    statistic <- ZI
    names(statistic) <- "Moran I statistic standard deviate"
    
    if (alternative == "two.sided") {
        PrI <- 2 * pnorm(abs(ZI), lower.tail = FALSE)
    } else if (alternative == "greater") {
        PrI <- pnorm(ZI, lower.tail = FALSE)
    } else {
        PrI <- pnorm(ZI)
    }
    
    if (!is.finite(PrI) || PrI < 0 || PrI > 1) {
        warning("Out-of-range p-value: reconsider test arguments")
    } 
    
    if (PrI == 0) {
        PrI <- scales::pvalue(PrI, accuracy = 0.00001)
    } else {
        PrI <- round(PrI, 6)
    }
    
    vec <- c(I, EI, VI)
    names(vec) <- c("Moran I statistic", "Expectation", "Variance")
    
    method <- paste("Moran I test under", ifelse(randomisation,
                                                 "randomisation", "normality"))
    
    data.name <- paste(xname, 
                       ifelse(rank, "using rank correction",""), 
                       "\nweights:", wname, 
                       ifelse(is.null(na.act), "", paste("\nomitted:", paste(na.act, collapse = ", "))),
                       ifelse(adjust.n && isTRUE(any(sum(spdep::card(listw$neighbours) == 0L))),
                              "n reduced by no-neighbour observations\n", ""),
                       ifelse(drop.EI2, "EI^2 term dropped in VI", ""), 
                       "\n")
    
    res <- list(statistic = statistic, p.value = PrI, estimate = vec, 
                alternative = alternative, method = method, data.name = data.name)
    
    if (!is.null(na.act)) {
        attr(res, "na.action") <- na.act
    }
    
    class(res) <- "htest"
    
    res
}

