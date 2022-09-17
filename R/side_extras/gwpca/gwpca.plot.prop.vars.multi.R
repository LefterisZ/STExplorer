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
        theme +
        theme(legend.position = "none")
}
