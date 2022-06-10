#' @description This function will gather the X/Y min and max from your set of 
#' points and will generate a bounding box (bbox). Then it will create and 
#' return a polygon (an sf object with geometry POLYGON).
#' 
#' The returned POLYGON is extended by one spot diameter to every direction. The
#' function will subtract from the xmin and ymin one spot diameter in pixels and
#' additionally will add to the xmax and ymax one spot diameter in pixels. This 
#' is done so that the polygon doesn't fall directly on the spots that are on 
#' the edges of the tissue.
#' 
#' @param input_points is an sf object that contains a set of points. The 
#' Geometry type must be of type POINT
#' @param sdiameter the diameter of each spot in pixels. Can be calculated by
#' using \code{sdiametereter()}
#' 
#' @export

make_bb_polygon <- function(input_points, sdiameter) {
    
    ## generate the bbox
    bbox <- st_bbox(input_points)
    
    ## get the bbox min/max dimensions into a 2x2 matrix
    points <- matrix(c(bb["xmin"]-sdiameter, bb["ymin"]-sdiameter, 
                       bb["xmin"]-sdiameter, bb["ymax"]+sdiameter,
                       bb["xmax"]+sdiameter, bb["ymax"]+sdiameter, 
                       bb["xmax"]+sdiameter, bb["ymin"]-sdiameter, 
                       bb["xmin"]-sdiameter, bb["ymin"]-sdiameter),
                     ncol = 2, byrow = T)
    
    ## Return
    st_polygon(list(points))
}
