#' Fetch gene expression for discrepancy location
#'
#' @name getDiscrepancyGeneData
#'
#' @description
#' A function to get the required data for the genes that make the spot in
#' question have a high discrepancy score. This function is based on the
#' gw.pcplot() function from the `GWmodel` package.
#'
#' @param m_sfe a SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' data object.
#'
#' @param assay Assay type for the spatial expression data (counts, logcounts,
#' etc.).
#'
#' @param vars Variables of interest (genes) to be evaluated. Default is NULL,
#' which includes all variables.
#'
#' @param focus Tissue locations of interest (barcodes) for which discrepancy
#' heatmaps will be generated.
#'
#' @param dMetric Distance metric used for generating the distance matrix which
#' will be used to identify the focus location's neighbors.
#'
#' @param sample_id Sample ID for which discrepancy data is retrieved.
#'
#' @param bw Bandwidth parameter for selecting neighbors for the heatmap.
#'
#' @param mean.diff Threshold for selecting genes based on the difference from
#' the mean discrepancy score. Default is 1.
#'
#' @param show.vars Display option for variables (genes) in the heatmap. Options
#' are "top" (genes with higher discrepancies) or "all" (all genes).
#'
#' @param exportExpression Defaults to FALSE. If TRUE, it prints out a data
#' frame of the gene data that were used in the discrepancy heatmap for the
#' specific discrepancy location. The data contain a column named 'gene_name'
#' which includes the gene IDs of the genes. ENSGene IDs are used as row names.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @importFrom dplyr mutate select filter across starts_with if_else
#' @importFrom dplyr left_join ends_with
#' @importFrom tibble column_to_rownames
#' @importFrom Matrix rowMeans
#' @importFrom SpatialFeatureExperiment colData rowData
#' @importFrom SummarizedExperiment assay
#'
#' @return
#' If exportExpression is FALSE, the function returns a list containing the
#' following:
#' - vars: Variables (genes) of interest.
#' - nbrlist: List of neighbours within the specified bandwidth.
#' - data.focus: Data frame containing the gene data for the discrep. heatmap.
#' - annot: Data frame containing the annotation data for the genes.
#' If exportExpression is TRUE, the function returns a data frame containing the
#' gene data that were used in the discrepancy heatmap.
#'
#' @examples
#' # Load data
#' data(sfe)
#'
#' # Set parameters
#' assay <- "logcounts"
#' focus <- colnames(sfe)[1] # outlier Barcode
#' sample_id <- "JBO019"
#' dMetric <- "euclidean"
#' bw <- 450
#' show.vars <- "top"
#'
#' # Get Gene expression data
#' discData <- getDiscrepancyGeneData(sfe, assay = assay, focus = focus,
#' dMetric = dMetric, sample_id = sample_id, bw = bw, show.vars = show.vars)
#' class(discData)
#' names(discData)
#'
#' @export
getDiscrepancyGeneData <- function(m_sfe,
                                   assay,
                                   vars = NULL,
                                   focus,
                                   dMetric,
                                   sample_id,
                                   diam,
                                   mean.diff = 1,
                                   show.vars = c("top", "all"),
                                   exportExpression = FALSE) {
    ## Check arguments
    show.vars <- match.arg(show.vars)

    ## SFE or metaSFE?
    sfe <- .int_sfeORmsfe(m_sfe = m_sfe)

    ## Get required data
    data <- assay(sfe, assay)
    annot <- rowData(sfe)[,c("id", "gene_name")]
    sample_locations <- colData(sfe)$sample_id == sample_id

    ## Set some variables
    dp.n <- ncol(data)  # Number of data points
    row.nm <- rownames(data)  # Gene IDs

    ## Check if variables of interest are provided
    if (is.null(vars)) {
        ## Use all variables if none are specified
        vars <- row.nm
    } else {
        ## Use the given variables but check that match the data
        var.check <- is.na(match(vars, row.nm))
        ## Check if vars input matches the data
        if (sum(var.check) > 0) {
          stop("`vars` input doesn't match with data.\n",
               "Please check that you provided a valid set of genes.")
        }
    }
    data <- data[vars, sample_locations]
    data <- as.matrix(data)
    annot <- as.data.frame(annot[vars,])

    ## Check if a distance matrix is given
    if (!missing(dMetric)) {
        dMat <- getDistMat(sfe = sfe, dMetric = dMetric)
        dMat <- dMat[sample_locations, sample_locations]
        ## Get the distances for the specific outlier in question
        dists <- dMat[,focus]
    } else {
      stop("A distance matrix is required.\n",
           "Please generate one using the addDistMat() function and provide ",
           "here the name of the metric used (i.e. `euclidean`)")
    }

    ## Fetch neighbours within the specified bandwidth
    nbrlist <- which(dists < diam + diam / 2)

    ## Prepare the data for the heatmap
    data.focus <- as.data.frame(data[,nbrlist])
    if (show.vars == "top") {
        data.focus <- data.focus %>%
            mutate(nb.mean = rowMeans(across(which(colnames(.) != focus))),
                   focus.diff = abs(.data$nb.mean - .data[[focus]])) %>%
            filter(.data$focus.diff > mean.diff) %>%
            select(-c("nb.mean", "focus.diff"))
    }

    data.focus <- data.focus %>%
        t() %>%
        as.data.frame() %>%
        mutate(Distance = round(dists[nbrlist], 0))

    ## Check that mean.diff cut off is not too strict.
    if (ncol(data.focus) == 1) {
      stop("Location with barcode: ", focus,
           "\nNo genes passed the cut off of differnce from mean set by the ",
           "`mean.diff` argument.\nYou can either:\n",
           "    1. lower the cut off.\n",
           "    2. if you provided a selection of genes in the `var` argument ",
           "\n       consider providing a different list of genes or leaving\n",
           "       the `var` argument as default (NULL) to include all genes.")
    }

    out <- list(vars = vars,
                nbrlist = nbrlist,
                data.focus = data.focus,
                annot = as.data.frame(annot))

    if (exportExpression) {
        out <- data.focus %>%
            dplyr::select(starts_with("ENSG")) %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column(var = "id") %>%
            left_join(annot) %>%
            mutate(gene_name = if_else(is.na(.data$gene_name),
                                       .data$id,
                                       .data$gene_name)) %>%
            column_to_rownames(var = "id")
    }

    return(out)
}

#' Get Discrepancy Data
#'
#' Retrieves discrepancy data based on sample ID from spatial expression data
#' and GWPCA results.
#' @name getDiscrepancyLocData
#'
#' @param m_sfe a SpatialFeatureExperiment or a MetaSpatialFeatureExperiment
#' data object.
#' @param gwpca GWPCA results object.
#' @param sample_id Sample ID for which discrepancy data is retrieved.
#'
#' @return A data frame containing the discrepancy data, including barcodes,
#' coordinates, discrepancy scores, and geometries.
#'
#' @details This function extracts the discrepancy data from spatial expression
#' data and GWPCA results based on a specified sample ID. It selects the sample
#' locations matching the given sample ID, identifies the discrepancies using
#' the is_outlier information from GWPCA, and retrieves the corresponding
#' barcodes, coordinates, discrepancy scores, and geometries. The discrepancy
#' data is returned as a data frame.
#'
#' @importFrom SpatialFeatureExperiment colData colGeometries spatialCoords
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Load data
#' data(sfe)
#' data(gwpca)
#'
#' # Set sample name
#' sample_id <- "JBO019"
#'
#' # Get Location data
#' discData <- getDiscrepancyLocData(sfe, gwpca, sample_id)
#' head(discData)
#'
#' @export
getDiscrepancyLocData <- function(m_sfe, gwpca, sample_id) {
  ## SFE or metaSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)
  #sample_locations <- colData(sfe)$sample_id == sample_id
  is_disc <- gwpca$CV$is_outlier

  discData <- data.frame(
      barcodes = colData(sfe)$Barcode[is_disc], #[sample_locations]
      coords = spatialCoords(sfe)[is_disc, ], # [sample_locations, ]
      discScore = gwpca$CV$CV[is_disc],
      geometry = colGeometry(sfe, "spotHex")[is_disc, ] # [sample_locations, ]
  )

  return(discData)
}



.int_sfeORmsfe <- function(m_sfe, sample_id) {
  SFE <- is(m_sfe, "SpatialFeatureExperiment")
  metaSFE <- is(m_sfe, "MetaSpatialFeatureExperiment")
  if (SFE) {
    sfe <- m_sfe
  } else if (metaSFE) {
    sfe <- getSFE(m_sfe, sample_id)
  }

  return(sfe)
}
