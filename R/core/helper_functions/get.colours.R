#' @name get.colours
#' 
#' @description 
#' 
#' @param 
#' 
#' @export

get.colours <- function(number){
    if (number > 167) {
        stop("The plot needs more than 167 colours. 
             Please consider plotting it using the default ggplot colour scheme.
             You can also provide a manual selection of more than 100 colours to
             the ggplot function.")
    }
    
    
    # Set some palettes in a list
    col.list <- list(palette36 = c4a("palette36"),
                     glasbey.32 = c4a("glasbey"),
                     alphabet2.26 = c4a("alphabet2"),
                     wright25 = c4a("wright25"),
                     light24 = c4a("light24"),
                     dark24 = c4a("dark24"))
    
    # Select the right number of colours needed
    if (number <= 25) {
        colours <- col.list$wright25[1:number]
    } else if (number <= 36) {
        colours <- col.list$palette36[1:number]
    } else if (number <= 50) {
        colours <- c(col.list$wright25,
                     col.list$glasbey.32[1:(number - 25)])
    } else if (number <= 75) {
        colours <- c(col.list$wright25,
                     col.list$light24,
                     col.list$alphabet2.26[1:(number - 49)])
    } else if (number <= 111) {
        colours <- c(col.list$wright25,
                     col.list$light24,
                     col.list$alphabet2.26,
                     col.list$palette36[1:(number - 75)])
    } else if (number <= 143) {
        colours <- c(col.list$wright25,
                     col.list$light24,
                     col.list$alphabet2.26,
                     col.list$palette36,
                     col.list$glasbey.32[1:(number - 111)])
    } else if (number <= 167) {
        colours <- c(col.list$wright25,
                     col.list$light24,
                     col.list$alphabet2.26,
                     col.list$palette36,
                     col.list$glasbey.32,
                     col.list$dark24[1:(number - 143)])
    }
    
    return(colours)
}