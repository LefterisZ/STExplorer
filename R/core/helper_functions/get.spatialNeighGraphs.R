#' @name get.spatialNeighGraphs
#' 
#' @description
#' A wraper around the \code{spdep} package functions that generate neighbour 
#' graphs and adds them in the \code{colGraphs} slot of the SFE 
#' (\code{SpatialFeatureExperiment}) object.
#' 
#' @param obj The SFE object.
#' 
#' @param sample_id character string specifying unique sample identifier for 
#' the sample you need to generate a neighbour graph.
#'  
#' @param type The different types that you can use to generate a neighbours 
#' graph. There are three categories; (a) Contiguity-based - "poly2nb", 
#' (b) Graph-based - "tri2nb", "soi.graph", "gabrielneigh", "relativeneigh" and
#' (c) Distance-based - "knearneigh", "dnearneigh". For more information about 
#' the individual functions visit the \code{spdep}'s documentation and vignette.
#' 
#' @param style default “raw”; style can take values “raw”, “W”, “B”, “C”, “U”, 
#' “minmax”, and “S” description. Argument passed to \code{spdep}'s 
#' \code{nb2listwdist} function. For more information about the individual 
#' functions visit the \code{spdep}'s documentation and vignette.
#' 
#' @param distMod default “idw”; the intended type of distance modelling, can 
#' take values “raw”, “idw”, “exp”, and “dpd”. Argument passed to \code{spdep}'s 
#' \code{nb2listwdist} function. For more information about the individual 
#' distance models visit the \code{spdep}'s documentation and vignette. The 
#' uses the \code{nb2listw} instead of the \code{nb2listwdist} function and it 
#' does not model the distance weights between the neighbours.
#' 
#' @param glist list of general weights corresponding to neighbours. Used only
#' when \code{(distMod == "raw")}
#' 
#' @param alpha default 1; a parameter for controlling the distance modelling. 
#' Argument passed to \code{spdep}'s \code{nb2listwdist} function, see 
#' \code{spdep}'s \code{nb2listwdist} function “Details” for more info.
#' 
#' @param dmax default NULL, maximum distance threshold that is required for 
#' weight type “dpd” but optional for all other types. Argument passed to 
#' \code{spdep}'s \code{nb2listwdist} function, see its help page for more info.
#' 
#' @param zero.policy	 default TRUE; if TRUE permit the weights list to be 
#' formed with zero-length weights vectors description, if FALSE stop with error 
#' for any empty neighbour sets, if NULL use global option value. Leave it as 
#' TRUE to avoid conflict with SFE object.
#' 
#' @param sym	 a logical argument indicating whether or not neighbours should be 
#' symmetric (if i->j then j->i) description. Used internally by the 
#' \code{graph2nb} function and only for graph-based neighbour types.
#' 
#' @param sfe a logical argument indicating whether or not the output should be 
#' added in the \code{colGraphs} slot of the SFE object or not. Default is TRUE;
#' it will be added; if FALSE a \code{listw} object is returned.
#' 
#' @param ... arguments that are passed down to the \code{spdep} functions 
#' called by the \code{type} argument. To see what else is needed for the functions to 
#' operate correctly visit the \code{spdep}'s documentation and vignette.
#' 
#' @export

get.spatialNeighGraphs <- function(obj,
                                   sample_id, 
                                   type = c("poly2nb", "tri2nb", "soi.graph",
                                            "gabrielneigh", "relativeneigh",
                                            "knearneigh", "dnearneigh"),
                                   style = c("raw", "W", "B", "C", "U", 
                                             "minmax", "S"),
                                   distMod = c("raw", "idw", "exp", "dpd"),
                                   glist = NULL,
                                   alpha = 1,
                                   dmax = NULL,
                                   zero.policy = TRUE,
                                   sym = FALSE,
                                   sfe = TRUE,
                                   ...) {
  ## Prepare some data
  type <- match.arg(type)
  style <- match.arg(style)
  distMod <- match.arg(distMod)
  data <- spatialCoords(obj)[colData(obj)$sample_id %in% sample_id, ]
  row.names <- colData(obj)[colData(obj)$sample_id %in% sample_id, "Barcode"]
  
  ## Generate the graph
  if (type == "poly2nb") {
    dataH <- colGeometries(obj)$spotHex[colData(obj)$sample_id %in% sample_id, ]
    nb_graph <- poly2nb(pl = dataH, row.names = row.names, ...)
    
  } else if (type == "tri2nb") {
    nb_graph <- tri2nb(coords = data, row.names = row.names)
    
  } else if (type == "soi.graph") {
    nb_tri <- tri2nb(coords = data, row.names = row.names)
    nb_graph <- graph2nb(soi.graph(tri.nb = nb_tri, coords = data, ...), 
                         row.names = row.names,
                         sym = FALSE)
    
  } else if (type == "gabrielneigh") {
    nb_graph <- graph2nb(gabrielneigh(coords = data, ...), 
                         row.names = row.names,
                         sym = FALSE)
    
  } else if (type == "relativeneigh") {
    nb_graph <- graph2nb(relativeneigh(coords = data, ...), 
                         row.names = row.names,
                         sym = FALSE)
    
  } else if (type == "knearneigh") {
    nb_graph <- knn2nb(knearneigh(x = data, ...), 
                       row.names = row.names,
                       sym = FALSE)
    
  } else if (type == "dnearneigh") {
    nb_graph <- dnearneigh(x = data, row.names = row.names, ...)
    
  }
  
  ## Make the simple graph a weighted list of neighbours 
  ## Check that a distance modelling has been selected.
  if (distMod == "raw") {
    nb_graph_w <- nb2listw(neighbours = nb_graph, glist = glist, 
                           style = style, zero.policy = zero.policy)
  } else {
    data <- colGeometries(obj)$spotCntd[colData(obj)$sample_id %in% sample_id, ]
    nb_graph_w <- nb2listwdist(neighbours = nb_graph, x = data, 
                               type = "idw", style = "W", 
                               zero.policy = zero.policy)
  }
  
  if (sfe) {
    colGraph(obj, sample_id) <- nb_graph_w
    obj
  } else if (!sfe) {
    nb_graph_w
  }
  
}