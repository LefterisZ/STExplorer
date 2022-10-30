#' @name gwpca.plot.outlier
#' 
#' @description A function to visualise the genes that make the spot in question
#'              have a high discrepancy score. This function is originating from
#'              the \code{gw.pcplot()} function of \code{GWmodel} package.
#' @param data  a data frame containing vst gene expression data.
#' @param vars  a vector of variable (gene) names to be evaluated. Default is 
#'              NULL and leading to evaluation of all variables (genes).
#' @param dMat  a distance matrix generated using \code{gw.dist()} function.
#' @param focus an integer pointing to the position of the location in focus 
#'              in the data input.
#' @param bw    bandwidth used in the weighting function.
#' @param mean.diff is the cut-off difference from the mean. To make the heatmap
#'                  more interpretable the function calculates for each gene the  
#'                  vst gene expression mean of the neighbours and from this  
#'                  subtracts the vst gene expression of the location in focus 
#'                  and takes the absolute values. By default, the cut-off value
#'                  is 1 meaning that a gene must have a difference > 1 to be 
#'                  present in the heatmap.
#' @param show.vars Takes two values "top" or "all". If "top" then the 
#'                  \code{mean.diff} argument is used to sort the genes shown in
#'                  the heatmap. If "all" then the \code{mean.diff} argument is
#'                  unused and all genes are shown in the heatmap. This is 
#'                  particularly useful when the user provides a small selection
#'                  of genes and doesn't want to further sort the list. Defaults
#'                  to "top".
#' @param scale pheatmap scale parameter. By default scale = "row" is used. For 
#'              more information look at the pheatmap documentation.
#' @param gene.names a TURE or FALSE value. Gives the option to print in the 
#'                   heatmap the gene names instead of the ENSGene IDs. Defaults 
#'                   to FALSE. If TRUE, it MUST be combined with a biomart table.
#' @param biomart a biomart table that includes a "gene_name" column from ENSembl
#'                BioMart annotation tables. It can be generated using the 
#'                \code{biomaRt} package or the \code{create_biomart} function.
#' @param show.data a TRUE or FALSE value. Returns a dataframe with the heatmap
#'                  data. Defaults to FALSE.
#' @param check.names a TRUE or FALSE value. It is passed to \code{data.frame}
#'                    inside \code{annotate_data}. Default is FALSE and if it is
#'                    TRUE it is only used when gene.names = TRUE as well. This
#'                    option is given because the check.names option can introduce
#'                    minor changes to the column names which then lead to errors
#'                    with the heatmap visualisation. A common change is "-" to 
#'                    ".". If check.names = FALSE then data.frame will not change
#'                    the column names. For more information look in 
#'                    \code{data.frame} help page.
#' @param ...   other graphical parameters, to be passed to \code{pheatmap} 
#'              function.
#' @export


gwpca.plot.outlier <- function(data, vars = NULL, dMat, focus, bw, 
                               mean.diff = 1, show.vars = "top", scale = "row", 
                               gene.names = FALSE, biomart = NULL,
                               check.names = FALSE, show.data = FALSE, ...) {

    #### Set some variables needed below ####
    dp.n <- nrow(data)              # Find the number of rows in data (data points)
    i <- focus                      # An integer pointing to the observation point (the outlier in question)
    col.nm <- colnames(data)        # Get the data column names (should be gene IDs)
    
    #### Check if user provided a certain set of genes to look at and prepare the data matrix ####
    if(is.null(vars)) {             # Check if variables of interest are provided
        vars <- colnames(data)      # IF NOT: use all variables
        var.idx <- 1:length(vars)   # Make the index for using all variables
        data <- as.matrix(data)     # Transform data to a matrix
    } else {                        # IF YES: use the given variables
        var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))] # Get the index of the variables to be evaluated
        data <- data[,var.idx]      # Select from data only the data for the variables to be evaluated
        data <- as.matrix(data)     # Transform it into a matrix
        if(length(var.idx) == 0){
            stop("Variables input doesn't match with data") # Check the length of this index. If == 0 then stop with error
        }
    }
    
    m <- ncol(data)                  # Get the number of columns (number of variables to be evaluated)
    
    #### Check if a distance matrix is given ####
    if (missing(dMat)) {             # Check if a distance matrix is given
        stop("A distance matrix is required. \n
             Please generate one using gw.dist() function")
    } else {                         # IF YES (else):
        dim.dMat<-dim(dMat)          # Get the dimensions of the given dMat
        if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n) { # Check if the dimensions match the dimensions of the dp.n
            stop("Dimensions of dMat are not correct") # IF NO: stop with an error
        }
        dists <- dMat[,i]            # Get the distances for the specific outlier in question
    }
    
    #### Find neighbours ####
    nbrlist <- which(dists < bw)  # Get the neighbours that are closer than the bandwidth

    #### Prepare the data for the heatmap ####
    focus.nm <- rownames(data)[i] # Get the name of the focus location
    if (show.vars == "top"){
        data.focus <- data[nbrlist,] %>% # get the expression data of the focus location and the neighbours
            t() %>%
            as.data.frame()  %>% 
            mutate(nb.mean = rowMeans(across(which(colnames(.) != focus.nm)))) %>% # find the vst mean of neighbours
            mutate(focus.diff = abs(.data$nb.mean - .data[[focus.nm]])) %>% # find absolute difference of mean to focus point 
            filter(.data$focus.diff > mean.diff) %>% # filter to the genes that have a difference from vst mean > mean.diff
            dplyr::select(-c("nb.mean", "focus.diff"))
    } else if (show.vars == "all") {
        data.focus <- data[nbrlist,] %>% # get the expression data of the focus location and the neighbours
            t() %>%
            as.data.frame()
    } else {
        stop("Please provide a valid value for show.vars. Either 'top' or 'all'")
    }
    
    #### Annotate the columns (tissue locations/ Barcodes) ####
    nbr.annotation <- data.frame(dists[nbrlist]) %>%
        mutate("Barcodes" = rownames(data)[nbrlist]) %>%
        column_to_rownames(var = "Barcodes") %>%
        rename("Distance" = colnames(.))

    #### Generate annotation colours based on distance from focus location ####
    nbr.lvls <- nbr.annotation %>%      # Find neighbour levels
        group_by(Distance) %>%
        tally()
    nbr.ann.colours <- scales::dscale(factor(1:count(nbr.lvls)[[1]]), 
                                      grey_pal(start = 0, end = 0.8))   # Get colour names from grey palette
    nbr.ann.colours <- rep(nbr.ann.colours, times = nbr.lvls$n)  # Update the colour names vector based on the levels
    nbr.ann.colours[1] <- "red"         # Make the focus point red. (since it has a distance of 0 it is always first)
    nbr.ann.colours <- list(Distance = nbr.ann.colours) # Convert it to a list
    
    #### Translate the ENSGene IDs to Gene names ####
    if(gene.names){
        if(is.null(biomart)){
            stop("Please provide a biomart table")
        } else {
            data.focus <- annotate_data.frame(data.focus, 
                                              biomart, 
                                              check.names = check.names)
            data.focus <- data.focus %>%
                column_to_rownames(var = "gene_name")
            data.focus2plot <- data.focus %>%
                dplyr::select(-c(1:3))
        }
    } else {
        data.focus2plot <- data.focus
    }
    
    #### Plot the Heatmap ####
    pheatmap(data.focus2plot,
             scale = scale,
             annotation_col = nbr.annotation,
             annotation_colors = nbr.ann.colours,
             ...)
    
    #### Provide the heatmap data ####
    if (show.data){
        return(data.focus)
    }
}


