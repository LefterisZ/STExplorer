#' @name gwpca.plot.prop.vars.multi
#' @description A function to create violin plots for Percentage of Total 
#'              Variation (PTV) from different number of Local PCs. can be used 
#'              alone and in an \code{lapply()} function to produce multiple plots.
#' @param gwpca.obj the gwpca object
#' @param n.comp the number of components to be used
#' @param violcol colour for the violin plots
#' 
#' @export

gwpca.plot.prop.vars.multi <- function(data, theme){
    # use melt from reshape2 (comes with R installation) to melt the data in a
    # two column df with variable and value colnames.
    ggplot(reshape2::melt(data), 
           aes(x = variable, y = value, fill = variable)) + 
        geom_violin(alpha = 0.8,
                    trim = FALSE) + 
        geom_boxplot(width = 0.1, 
                     fill = "grey40", 
                     outlier.colour = NA) +
        ylab("Percentage of Total Variation (PTV)") +
        xlab("") +
        ylim(0, 100) +
        theme +
        theme(legend.position = "none")
}

#' @name gwpca.plot.prop.vars.single
#' @description A function to create violin plots for Percentage of Total 
#'              Variation (PTV) from different number of Local PCs. can be used 
#'              alone and in an \code{lapply()} function to produce multiple 
#'              plots.
#' @param gwpca.obj the gwpca object
#' @param n.comp the number of components to be used
#' @param violcol colour for the violin plots
#' 
#' 
#' 
#' @export

gwpca.plot.prop.vars.single <- function(data, column = NULL, violcol = "skyblue", 
                                        dotcol = "white", boxcol = "grey40", 
                                        theme = my_theme, limY = NULL, 
                                        dotsize = 0.7, xlab = NULL, 
                                        ylab = "Percentage of Total Variation (PTV)"){
    # use melt from reshape2 (comes with R installation) to melt the data in a
    # two column df with variable and value colnames.
    # CONSIDER pivot_longer from tidyverse instead of melt
    
    # Generate the internal plot function in.plot() as the main plot function to 
    #   avoid repeating the code in the if-else statements.
    in.plot <- function(dat, violcol, dotcol, boxcol, theme, limY, xlab, ylab,
                        dotsize){
        # Check if certain Y-axis limits are set
        if (is.null(limY)) {
            limY <- c(min(dat$value) - 5, max(dat$value) + 5)
        }
        
        # Plot
        ggplot(dat, 
               aes(x = variable, y = value, fill = variable)) + 
            geom_violin(alpha = 0.8,
                        fill = violcol,
                        trim = FALSE) + 
            ggbeeswarm::geom_beeswarm(colour = dotcol,
                                      size = dotsize,
                                      cex = 1) +
            geom_boxplot(width = 0.1, 
                         fill = boxcol, 
                         outlier.colour = NA) +
            ylab(ylab) +
            xlab(xlab) +
            ylim(limY) +
            theme +
            theme(legend.position = "none")
    }
    
    
    # Run this if else statement to check if the function is used inside an 
    #   apply() function.
    if (is.null(column)){
        dt <- reshape2::melt(data)
        in.plot(dat = dt, violcol, dotcol, boxcol, theme, limY, xlab, ylab,
                dotsize)
    } else {
        dt <- reshape2::melt(data[, column])
        dt <- mutate(dt, variable = column)
        in.plot(dat = dt, violcol, dotcol, boxcol, theme, limY, xlab, ylab,
                dotsize)
    }
}

