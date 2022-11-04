#' @name poly2nb.ste
#' 
#' @description a function to find nbs using a polygons tessellation. It
#'              builds arround \code{poly2nb} function from spdep package. It 
#'              adds the ability to identify polygons that have no nbs 
#'              and remove them from the final list.
#'
#' @param pls list of polygons of class extending SpatialPolygons, or an sf or 
#'            sfc object containing non-empty (multi-)polygon objects.
#' @param snap boundary points less than snap distance apart are considered to 
#'             indicate contiguity.
#' @param remove a TRUE or FALSE value. Defaults to TRUE. The function will 
#'               remove all spots with no neighbours.
#' @param show a TRUE or FALSE value. Defaults to FALSE The function will print 
#'             the XY coordinates of the neighbourless locations. This requires
#'             the object supplied to the pls parameter to include the Barcode, 
#'             pixel_x and pixel_y columns apart from the polygon geometries.
#' @param columns a vector containing 3 columns. One column must be the region ID
#'                (Barcode by default), one must be the X coordinates 
#'                (pixel_x by default) and one must be the Y coordinates
#'                (pixel_y by default).
#' @param ... extra arguments to be passed into poly2nb. See ?spdep::poly2nb 
#'            documentation for more information.
#' 
#' @export


poly2nb.ste <- function(pls, snap = 0, remove = TRUE, show = FALSE, 
                        columns = c("Barcode", "pixel_x", "pixel_y"), ...){
    
    # Find the neighbours
    nbs <- poly2nb(pls, snap, ...)
    
    # Add region IDs (Barcodes) to the sub-lists
    names(nbs) <- attr(nbs, "region.id")
    
    # Find the locations with no neighbours
    no.nb.1 <- names(nbs)[unlist(lapply(nbs, function(x) sum(x) == 0))]
    if(is_empty(no.nb.1)) {
        message("Every location in the dataset has at least 1 neighbour")
    } else if(!is_empty(no.nb.1)) {
        message(length(no.nb.1), " locations in the dataset have no neighbours:")
        print(no.nb.1)
    }
    
    # Show the coordinates of the neighbourless locations
    if(show){
        check <- sum(columns %in% colnames(pls))
        if(check == 3){
            no.nb.2 <- filter(pls, Barcode %in% no.nb.1) %>%
                st_drop_geometry() %>%
                select(all_of(columns))
            print(no.nb.2)
        } else {
            stop("You've set 'show = TRUE' but the object provided in the 'pls' 
                 parameter does not contain all columns needed. Have a look at
                 ?poly2nb.ste for more information. If you provided your own 
                 column names then make sure there are no typos or missing 
                 columns.")
        }
    }
    
    # Remove the neighbourless locations
    if(remove){
        
        # Remove the neighbourless locations from the pls input
        pls.2 <- pls[!rownames(pls) %in% no.nb.1,]
        
        # Recalculate the neighbours
        nbs <- poly2nb(pls.2, snap, ...)
        
        # Set the extra required attributes
        attr(nbs, "no.nbs") <- no.nb.1 # set the neighbourless locations
    }
    
    # Return neighbours object
    return(nbs)
} 


