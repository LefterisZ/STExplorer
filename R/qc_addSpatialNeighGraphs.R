#' Add Neighbour Graphs
#'
#' @name addSpatialNeighGraphs
#'
#' @description
#' A function to add spatial neighbour graphs to a SpatialFeaturesExperiment
#' object. The function generates different types of spatial neighbour graphs
#' based on the coordinates of the spatial features (spots) in the object. The
#' generated graphs can be either raw neighbour graphs or weighted neighbour
#' graphs based on specified distance models. The function supports various
#' types of neighbour graph algorithms and distance models.
#' The function wraps around `spdep`'s neighbour functions.
#'
#' @param msfe The metaSFE object.
#'
#' @param sample_id character string, TRUE or NULL specifying sample/image
#' identifier(s); here, TRUE is equivalent to all samples/images.
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
#' distance models visit the \code{spdep}'s documentation and vignette. It
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
#' formed with zero-length weights vectors description, if FALSE stop with
#' error for any empty neighbour sets, if NULL use global option value. Leave
#' it as TRUE to avoid conflict with SFE object.
#'
#' @param sym	 a logical argument indicating whether or not neighbours should
#' be symmetric (if i->j then j->i) description. Used internally by the
#' \code{graph2nb} function and only for graph-based neighbour types.
#'
#' @param sfe_out a logical argument indicating whether or not the output
#' should be added in the \code{colGraphs} slot of the SFE object or not.
#' Default is TRUE; it will be added; if FALSE a \code{listw} object is
#' returned.
#'
#' @param ... arguments that are passed down to the \code{spdep} functions
#' called by the \code{type} argument. To see what else is needed for the
#' functions to operate correctly visit the \code{spdep}'s documentation and
#' vignette.
#'
#' @details This function adds spatial neighbour graphs to a
#' SpatialFeatureExperiment object. The neighbour graphs provide information
#' about the spatial relationships between the features (genes) in the object.
#' The neighbour graphs can be either raw neighbour graphs or weighted neighbour
#' graphs. The weighted neighbour graphs assign weights to the edges based on
#' distances or other measures, providing a more detailed representation of the
#' spatial relationships.
#'
#' @importFrom spdep poly2nb tri2nb soi.graph gabrielneigh relativeneigh
#' @importFrom spdep nb2listw nb2listwdist knearneigh dnearneigh knn2nb
#' @importFrom spdep graph2nb
#' @importFrom SpatialFeatureExperiment colData colGeometries spatialCoords
#' @importFrom SpatialFeatureExperiment colGraph colGraph<-
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @returns
#' If \code{sfe} is set to \code{TRUE}, the function returns the updated
#' SpatialFeatureExperiment object with the added neighbour graph stored in
#' the colGraphs slot. If \code{sfe} is set to \code{FALSE}, the function
#' returns the generated neighbour graph as a weighted list of neighbours.
#'
#' @examples
#' # Load a SpatialFeatureExperiment object
#' data(sfe)
#'
#' # Add spatial neighbour graphs in the provided SFE object
#' sfe <- addSpatialNeighGraphs(sfe = sfe, sample_id = "JBO019",
#' type = "knearneigh", style = "W", distMod = "raw", k = 6, sfe_out = TRUE)
#'
#' # Get the generated neighbour graph as an object of its own
#' nbr_graph <- addSpatialNeighGraphs(sfe = sfe, sample_id = "JBO019",
#' type = "knearneigh", style = "W", distMod = "raw", k = 6, sfe_out = FALSE)
#'
#' @export
addSpatialNeighGraphs <- function(msfe,
                                  sample_id = TRUE,
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
                                  sfe_out = TRUE,
                                  ...) {
  ## check arguments
  # stopifnot(is(msfe, "SpatialFeatureExperiment"))
  type <- match.arg(type)
  style <- match.arg(style)
  distMod <- match.arg(distMod)

  ## Select samples
  ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)

  ## Generate the graphs
  msfe_int <- lapply(msfe@sfe_data[ids], .int_addSpNghGphs,
                     # sample_id,
                     type = type,
                     style = style,
                     distMod = distMod,
                     glist = glist,
                     alpha = alpha,
                     dmax = dmax,
                     zero.policy = zero.policy,
                     sym = sym,
                     sfe_out = sfe_out,
                     ...)

  ## If specific samples where modified replace in the metaSFE list
  if (is.character(sample_id)) {
    msfe@sfe_data[names(msfe_int)] <- msfe_int
  } else {
    msfe@sfe_data <- msfe_int
  }

  return(msfe)
}



.int_addSpNghGphs <- function(sfe,
                              # sample_id,
                              type, #= c("poly2nb", "tri2nb", "soi.graph",
                                   #   "gabrielneigh", "relativeneigh",
                                   #   "knearneigh", "dnearneigh"),
                              style,# = c("raw", "W", "B", "C", "U",
                                     #  "minmax", "S"),
                              distMod,# = c("raw", "idw", "exp", "dpd"),
                              glist,# = NULL,
                              alpha,# = 1,
                              dmax,# = NULL,
                              zero.policy,# = TRUE,
                              sym,# = FALSE,
                              sfe_out,# = TRUE,
                              ...) {
  ## Prepare some data
  # type <- match.arg(type)
  # style <- match.arg(style)
  # distMod <- match.arg(distMod)
  data <- spatialCoords(sfe)#[colData(sfe)$sample_id %in% sample_id, ]
  row.names <- colData(sfe)[, "Barcode"]#colData(sfe)$sample_id %in% sample_id

  ## Generate the graph
  if (type == "poly2nb") {
    dataH <- colGeometries(sfe)$spotHex#[colData(sfe)$sample_id %in% sample_id,]
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
    data <- colGeometries(sfe)$spotCntd#[colData(sfe)$sample_id %in% sample_id, ]
    nb_graph_w <- nb2listwdist(neighbours = nb_graph, x = data,
                               type = "idw", style = "W",
                               zero.policy = zero.policy)
  }

  if (sfe_out) {
    colGraph(sfe) <- nb_graph_w #, sample_id
    sfe
  } else if (!sfe_out) {
    nb_graph_w
  }

}
