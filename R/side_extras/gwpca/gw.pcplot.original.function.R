gw.pcplot <- function(data,vars,focus,bw,adaptive = FALSE, ylim=NULL,ylab="",fixtrans=FALSE, p=2, theta=0, longlat=F,dMat,...) 
{
    if (is(data, "Spatial"))      # expects a spatial df with data and coordinates
    {                             # IF TRUE:
        p4s <- proj4string(data)  # Get data to p4s 
        loc<-coordinates(data)    # Get coordinates to loc
    }
    else                          # IF FLASE:
        stop("Given data must be a Spatial*DataFrame") # Stop with an error
    data <- as(data, "data.frame")# Make data a data.frame
    dp.n<-nrow(data)              # Find the number of rows in data (data points)
    i<-focus                      # An integer pointing to the observation point (the outlier in question)
    col.nm<-colnames(data)        # Get the data column names (should be gene IDs)
    var.idx<-match(vars, col.nm)[!is.na(match(vars, col.nm))]              # Get the index of the variables to be evaluated from the column names object from above
    if(length(var.idx)==0) stop("Variables input doesn't match with data") # Check the length of this index. If == 0 then stop with error
    x<-data[,var.idx]             # Select from data only the data for the variables to be evaluated
    x<-as.matrix(x)               # Transform it into a matrix
    m <- ncol(x)                  # Get the number of columns (number of variables to be evaluated)
    if (missing(dMat))            # Check if a distance matrix is given
    {                             # IF NOT:
        DM.given<-F                     # Set DM.given to FALSE
        if(dp.n <= 5000)                # Check if the data points are less than 5000
        {                               # IF YES:
            dMat <- gw.dist(dp.locat=loc, p=2, theta=0, longlat=FALSE) # Use the gw.dist function to generate the distance matrix
            DM.given<-T                 # Set DM.given to TRUE
        }
    }
    else                          # IF YES (else):
    {
        DM.given<-T                    # Set DM.given to TRUE
        dim.dMat<-dim(dMat)            # Get the dimensions of the given dMat
        if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n) # Check if the dimensions match the dimensions of the dp.n
            stop("Dimensions of dMat are not correct") # IF NO: stop with an error
    }
    
    if(DM.given)                  # Check if DM.given is TRUE. IF IT IS:
        dists <- dMat[,i]         # Get the distances for the specific outlier in question
    else                          # IF IT IS NOT:
        dists <-gw.dist(dp.locat=loc, focus = i, p=p, theta=theta, longlat=longlat) # Use gw.dist function to generate the distance matrix for the specific outlier in question
    if(adaptive)                  # IF IT IS SET TO ADAPTIVE: Calculate the distances
    {
        rnk<- rank(dists, ties.method='first')
        bw<- dists[rnk==bw]  
    }
    #dists <- dists**2
    nbrlist <- which(dists < bw)  # Get the neighbours that are closer than the bandwidth
    dists <- dists**2             # Raise the distances to the power of 2
    wts <- (1 - dists/(bw*bw))^12 # Calculate the weights for the distances
    #  if(adaptive)
    #    bw<-dists[bw]
    #nbrlist <- which(dists < bw*bw)
    xss <- scale(x)               # Scale the data to be easier ploted
    span <- 1:m                   # Get the span from 1 to the number of genes to be evaluated
    tsc <- 25/length(nbrlist)     # Calculate a factor that will give the different colouring on the plot based on the distance from the data point in question
    if (is.null(ylim))            # Plot
        ylim <- c(min(xss[nbrlist,]),max(xss[nbrlist,]))
    plot(span,xss[i,],type='l',ylim=ylim,
         xlim=c(0.5,m+0.5),col='red',lwd=6,axes=FALSE,xlab="",ylab=ylab) # put back after ylab: ,...
    #axis(1,at=1:m,labels=colnames(x),las=2,cex.axis=1.2)
    axis(2,at=seq(floor(ylim[1]),ceiling(ylim[2]),by=1),cex.axis=1.2)
    #abline(v=1:m,col=grey(0.6))
    lines(c(1,m),c(0,0),col=grey(0.6))
    if (fixtrans) {
        for (nbr in nbrlist) 
            lines(span,xss[nbr,],col=rgb(0,0,0,0.3), lwd=3) }
    else 
    {
        for (nbr in nbrlist)
            lines(span,xss[nbr,],col=rgb(0,0,0,tsc*wts[nbr]), lwd=3)  
    } 
    
}


data <- as.data.frame(data.mat) # vst data (gene expression data)
bw <- 3*spot_diameter(spatialDir)
i <- which(discrepancy_df == max(discrepancy_df))
st_geometry(polygons) <- "geom_cntd"
loc <- st_coordinates(polygons)
dMat <- dist.Mat

focus.nm <- rownames(data)[i]

data.focus <- data[nbrlist,] %>% # get the expression data of the focus location and the neighbours
    t() %>%
    as.data.frame() %>%
    mutate(nb.mean = rowMeans(across(which(colnames(.) != focus.nm)))) %>% # find the vst mean of neighbours
    mutate(focus.diff = abs(.data$nb.mean - .data[[focus.nm]])) %>% # find absolute difference of mean to focus point 
    mutate(labels = ifelse(.data$focus.diff > 1, TRUE, FALSE)) # add labels to the genes that have a difference from vst mean > 1


nb.alphas <- tsc*wts[nbrlist[!nbrlist %in% i]] # get the alpha values for the nbrs only

# create the input data table
inputPCAgw.ex.outlier <- counts[select,] %>% # select top 500 variable genes
    t() %>%                                  # transpose
    as.data.frame() %>%                      # make it a df
    rownames_to_column(var = "Barcode") %>%  # get barcodes in a column
    merge(., polygons[,c("Barcode")]) %>% # merge with coordinates
    column_to_rownames(var = "Barcode") %>%  # return barcodes to row names
    .[nb_names,]                             # order rows
st_geometry(inputPCAgw.ex.outlier) <- "geom_cntd" # make geom_cntd the geometry column
    

# Plot extreme outlier with labels for genes with difference from the mean of more than 1
xss.focus <- data.frame(focus.loc = xss[i,]) # get the scaled expression of the location in focus
xss.focus.nbrs <- xss[nbrlist[!nbrlist %in% i],] %>% # get the scaled expression data for the neighbours of the focus location
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_IDs") %>%
    reshape2::melt(id = "gene_IDs")

ggplot() +
    geom_line(data = xss.focus, aes(x = span, y = focus.loc), 
               size = 1.5, colour = "red") + 
    geom_line(data = xss.focus.nbrs,
              aes(x = gene_IDs,
                  y = value,
                  colour = variable))
    


