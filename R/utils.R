#' Add Distance Matrix
#'
#' Calculate and store the distance matrix based on spatial coordinates using
#' the specified metric.
#' @name addDistMat
#'
#' @param msfe \code{MetaSpatialFeatureExperiment} (Meta-SFE) object.
#' @param p Numeric scalar specifying the power parameter for the Minkowski
#' distance. See details below for more information.
#' @param sample_id character string or TRUE specifying sample/image
#' identifier(s); here, TRUE is equivalent to all samples/images.
#' @param ... Additional arguments to be passed to \code{gw.dist} function.
#'
#' @importFrom GWmodel gw.dist
#'
#' @return An updated SpatialFeatureExperiment object with the distance matrix
#' stored in the metadata.
#'
#' @details
#' This function calculates the distance matrix based on the spatial coordinates
#' of the provided SpatialFeatureExperiment object. The distance matrix is
#' computed using the specified metric, where p = 2 corresponds to Euclidean
#' distance, p = 1 corresponds to Manhattan distance, and p > 2 corresponds to
#' Minkowski distance with the specified power. The resulting distance matrix is
#' stored in the metadata of the SpatialFeatureExperiment object.
#'
#' @seealso [getDistMat()] for a function to get the distance matrix from within
#' an SFE object.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Load sfe object
#' data(sfe)
#'
#' # Add distance matrix to sfe
#' sfe <- addDistMat(sfe, p = 2)
#'
#'
#' @export
addDistMat <- function(msfe, p, sample_id = TRUE, ...) {
  ## Check sfe or msfe
  if (is(msfe, "SpatialFeatureExperiment")) {
    sfe <- TRUE
  } else if (is(msfe, "MetaSpatialFeatureExperiment")) {
    sfe <- FALSE
  }

  if (sfe) {
    sfe <- .int_addDistMat(sfe = msfe, p = p, ...)
    return(sfe)

  } else {
    ## Select samples
    ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)
    ## Generate the distance matrices
    msfe_int <- lapply(msfe@sfe_data[ids], .int_addDistMat, p = p, ...)
    ## If specific samples where modified replace in the metaSFE list
    if (is.character(sample_id)) {
      msfe@sfe_data[names(msfe_int)] <- msfe_int
    } else {
      msfe@sfe_data <- msfe_int
    }
    return(msfe)
  }
}

#' Get Distance Matrix
#'
#' Retrieve the distance matrix based on the specified metric from the metadata
#' of a SpatialFeatureExperiment object.
#'
#' @param msfe \code{MetaSpatialFeatureExperiment} (Meta-SFE) object.
#' @param dMetric Character string specifying the distance metric.
#' @param sample_id character string, or TRUE specifying sample/image
#' identifier(s); here, TRUE is equivalent to all samples/images.
#'
#' @return Distance matrix corresponding to the specified metric.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Load sfe object
#' data(sfe)
#'
#' # Add distance matrix to sfe
#' sfe <- addDistMat(sfe, p = 2)
#'
#' # Fetch distance matrix from sfe
#' dMat <- getDistMat(sfe, dMetric = "euclidean")
#'
#' # Check distance matrix
#' dMat[1:4, 1:4]
#'
#' @export
getDistMat <- function(msfe, dMetric, sample_id = TRUE) {
  ## Check sfe or msfe
  if (is(msfe, "SpatialFeatureExperiment")) {
    sfe <- TRUE
  } else if (is(msfe, "MetaSpatialFeatureExperiment")) {
    sfe <- FALSE
  }

  if (sfe) {
    dMat_int <- .int_getDistMat(msfe, dMetric = dMetric)
  } else {
    ## Select samples
    ids <- .int_getMSFEsmplID(msfe = msfe, sample_id = sample_id)
    dMat_int <- lapply(msfe@sfe_data[ids], .int_getDistMat, dMetric = dMetric)
  }

  return(dMat_int)
}


#' Get a vector of colours for plotting
#'
#' This function returns a vector of colours from pre-defined palettes based on
#' the desired number of colours.
#'
#' @param number The desired number of colours.
#'
#' @return A vector of colours.
#'
#' @details This function selects colours from various pre-defined palettes to
#' accommodate the desired number of colours. It supports up to 167 colours. If
#' a number greater than 167 is provided, an error is thrown.
#'
#' @importFrom cols4all c4a
#'
#' @seealso [c4a()] from the \code{\link{cols4all}} package.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Get a set of 10 colours
#' colours <- getColours(10)
#' colours
#'
#' @export

getColours <- function(number){
  if (number > 167) {
    stop("The plot needs more than 167 colours.\n",
         "Please consider plotting it using the default ggplot ",
         "colour scheme.\n",
         "You can also provide a manual selection of more than ",
         "100 colours to the ggplot function.")
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

#' Get outlier cut off value
#'
#' @name outlierCutoff
#' @description
#' This function identifies outliers in a vector of data based on a coefficient
#' and the interquartile range (IQR).
#'
#' @param dt A numeric vector of data.
#' @param coef The coefficient used to determine the cutoff for outliers.
#' The default value is 1.5, which is commonly used in boxplots to define
#' outliers as values that are more than 1.5 times the IQR away from the upper
#' or lower quartiles.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{out_up}: The upper cutoff value for outliers.
#'   \item \code{out_down}: The lower cutoff value for outliers.
#'   \item \code{outliers}: A logical vector indicating whether each data point
#'   is an outlier or not. \code{TRUE} indicates an outlier, and \code{FALSE}
#'   indicates a non-outlier.
#' }
#'
#' @importFrom stats quantile
#'
#' @examples
#' # Generate a vector of random data
#' data <- rnorm(100)
#'
#' # Detect outliers using the default coefficient (1.5)
#' outlierCutoff(data)
#'
#' # Detect outliers using a custom coefficient (2.0)
#' outlierCutoff(data, coef = 2.0)
#'
#' @export
outlierCutoff <- function(dt, coef = 1.5) {
  ## Calculate the quartiles of the data
  quantiles <- quantile(dt, probs = c(0.25, 0.75))

  ## Calculate the interquartile range (IQR)
  IQR <- quantiles[2] - quantiles[1]

  ## Calculate the cutoff values for outliers
  out_up <- quantiles[2] + coef * IQR
  out_down <- quantiles[1] - coef * IQR

  ## Identify outliers based on the cutoff values
  outliers <- dt < out_down | dt > out_up

  ## Create a list with the results
  res <- list(out_up = out_up, out_down = out_down, outliers = outliers)

  return(res)
}


#' Generate a Gradient of Colours
#'
#' This function creates a gradient of colours from a given vector of colours to white,
#' with a specified number of steps.
#'
#' @param colours A character vector representing the colours to use for the gradient.
#' @param steps An integer specifying the number of steps in the gradient.
#'
#' @return A list of vectors, each containing the colour values representing the gradient for a single input colour.
#'
#' @rdname getGradients
#'
#' @examples
#' getGradients(c("red", "green", "blue"), steps = 3)
#' getGradients(c("orange", "purple"), steps = 5)
#'
#' @export
getGradients <- function(colours, steps) {
  # Validate input parameters
  if (!is.character(colours)) {
    stop("Colours must be a character vector")
  }
  if (!is.numeric(steps) || steps <= 0) {
    stop("Steps must be a positive integer")
  }

  # Create a gradient for each input colour
  gradients <- lapply(colours, .int_getGradient, steps = steps)
  return(unlist(gradients))
}


#' Read data from a Curio Seeker (Slide-seq) experiment
#'
#' Reads spatial transcriptomics data from specified directory paths and
#' formats it into a `SpatialFeatureExperiment` object. This function supports
#' data in AnnData, Seurat, or CSV formats.
#'
#' @param samples A character vector specifying one or more directories, each
#'        corresponding to a 10x Genomics Visium sample. These directories
#'        should contain the count data files in one of the supported formats
#'        (AnnData, Seurat, CSV). If the `samples` parameter is named, those
#'        names will be used as sample identifiers.
#' @param sample_id A character vector specifying unique sample identifiers,
#'        one for each directory specified via `samples`. This parameter is
#'        ignored if `samples` is a named vector.
#' @param fileType A character string specifying the type of format to read
#'        count data from. Possible values are "AnnData", "Seurat", and "CSV".
#'        Defaults to the first element of the vector ("AnnData", "Seurat",
#'        "CSV"). More info in the details section below.
#'
#' @details
#' The function processes directories specified in the `samples` parameter, each
#' expected to represent a sample in the spatial transcriptomics experiment.
#' Based on the `fileType` parameter, the function looks for specific data
#' files:
#' \itemize{
#'   \item \strong{AnnData}: Looks for `.h5ad` files containing the count
#'         matrix, feature data, and spatial coordinates.
#'   \item \strong{Seurat}: Expects `.rds` files which are R serialised objects
#'         of Seurat objects containing similar data as above.
#'   \item \strong{CSV}: Reads raw CSV files for counts, and TSV
#'         files for barcodes and features, along with a CSV for spatial
#'         coordinates. The counts csv MUST include `counts.csv` in its name.
#' }
#'
#' For `AnnData` format, the function  uses the `anndata` package to load the
#' data.
#'
#' For `Seurat` format, it expects a Seurat object with counts in
#' `obj@assays$RNA@counts`, and spatial coordinates in
#' `obj@images$slice1@coordinates`.
#'
#' For `CSV` format, it expects files that contain a `counts.csv` in their name
#' for the counts matrix, `barcodes.tsv` for barcode names, `genes.tsv` for
#' feature names, and `MatchedBeadLocation.csv` for spatial coordinates.
#'
#' @return A `SpatialFeatureExperiment` object containing the spatial
#'         transcriptomics data with counts, spatial coordinates, and feature
#'         data.
#'
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom utils read.csv
#' @importFrom SpatialFeatureExperiment SpatialFeatureExperiment
#' @importFrom anndata read_h5ad
#' @importFrom Matrix as.matrix
#'
#' @examples
#' \dontrun{
#'   sfe <- readCurioSeeker(samples = c("sample1/", "sample2/"),
#'                          sample_id = c("sample1", "sample2"),
#'                          fileType = "AnnData")
#' }
#'
#' @export
#' @param samples a character vector specifying one or more directories, each
#' corresponding to a 10x Genomics Visium sample (see Details); if provided,
#' names will be used as sample identifiers
#' @param sample_id character string specifying unique sample identifiers, one
#' for each directory specified via samples; ignored if !is.null(names(samples))
#' @param fileType character string specifying the type of format to read count
#' data from. One of (AnnData, Seurat, CSV).
#'
readCurioSeeker <- function(samples,
                            sample_id,
                            fileType = c("AnnData", "Seurat", "CSV")) {
  if (fileType == "AnnData") {
    filePath <- .int_findFileWithRegex(directory = samples,
                                       pattern = "anndata|h5ad")
    obj <- anndata::read_h5ad(filename = filePath)

    ## Format counts table
    counts <- obj$X %>% # Matrix::as.matrix(sparse = TRUE) %>%
      Matrix::t() # %>%
      # as(Class = "dgCMatrix")

    ## Load coordinates
    coords <- obj$obsm$X_spatial %>%
      as.data.frame() %>%
      dplyr::mutate("Barcode" = rownames(obj))
    colnames(coords) <- c("x", "y", "Barcode")

    ## Create rowData
    rowData <- data.frame(symbol = colnames(obj))

  } else if (fileType == "Seurat") {
    filePath <- .int_findFileWithRegex(directory = samples,
                                       pattern = "seurat|rds")
    obj <- readRDS(filePath)

    ## Format counts table
    counts <- obj@assays$RNA@counts

    ## Load coordinates
    coords <- obj@images$slice1@coordinates %>%
      rownames_to_column(var = "Barcode")

    ## Create rowData
    rowData <- data.frame(symbol = rownames(obj))

  } else if (fileType == "CSV") {
    filePath <- .int_findFileWithRegex(directory = samples,
                                       pattern = "counts.csv")
    filePath_barcodes <- .int_findFileWithRegex(directory = samples,
                                                pattern = "barcodes.tsv")
    filePath_features <- .int_findFileWithRegex(directory = samples,
                                                pattern = "genes.tsv")
    filePath_coords <- .int_findFileWithRegex(directory = samples,
                                              pattern = "MatchedBeadLocation.csv")

    ## Load and format counts table
    barcodes <- read.table(filePath_barcodes)[,1]
    features <- read.table(filePath_features)[,1]
    counts <- read.csv(filePath, header = TRUE)
    count_dims <- dim(counts)
    if (count_dims[1] == length(barcodes)) {
      counts <- t(counts)
    }
    rownames(counts) <- features
    colnames(counts) <- barcodes

    ## Load coordinates
    coords <- read.csv(filePath_coords, header = TRUE)
    colnames(coords) <- c("Barcode", "x", "y")

    ## Create rowData
    rowData <- data.frame(symbol = features)
  }

  sfe_out <- SpatialFeatureExperiment(list(counts = counts),
                                      colData = coords,
                                      rowData = rowData,
                                      sample_id = sample_id,
                                      spatialCoordsNames = c("x", "y"),
                                      spotDiameter = NA_real_,
                                      annotGeometryType = "POINT")

  return(sfe_out)
}

#' Load and Process Visium SpatialFeatureExperiment (SFE) Data
#'
#' This function loads and processes SpatialFeatureExperiment (SFE) data from
#' within the package's `extdata` folder. It addresses the issue of missing
#' image data when using the `usethis::use_data` function to store SFE objects
#' in `.rda` format. The function ensures that the data is loaded correctly and
#' that the images are reloaded from the external data located in the package's
#' `inst/extdata/visium/` directory.
#'
#' Depending on the dataset specified, the function loads relevant data objects
#' (`msfe_2` for lung data or `msfe` for prostate data), processes each sample
#' by adding the corresponding low-resolution spatial tissue images, and applies
#' the appropriate scale factors. It also mirrors the images as part of the
#' processing.
#'
#' This function is essential for working with spatial data, as it reconstructs
#' the image pointers (which are lost when saving SFE objects via
#' `usethis::use_data`) and ensures the correct association between images and
#' sample data.
#'
#' @param dataset A character string specifying the dataset to load. Must be
#'   either `"lung"` or `"prostate"`. Based on this input, the function loads
#'   the corresponding SFE object and processes the associated data.
#'
#' @details
#' For the `"lung"` dataset, the function loads the `msfe_2` object and
#' processes samples located in the `inst/extdata/visium/lung/Fibrotic/` and
#' `inst/extdata/visium/lung/Healthy/` directories. For the `"prostate"`
#' dataset, it loads the `msfe` object and processes the sample located in
#' `inst/extdata/visium/prostate/H2_5/`.
#'
#' The function retrieves scale factors from the `scalefactors_json.json` file
#' within each sample's directory and adds the corresponding low-resolution
#' tissue image to the SFE object. It also mirrors the image using the
#' `SpatialFeatureExperiment::mirrorImg` function.
#'
#' @return Returns the updated SFE object after processing, with images properly
#' reloaded and associated with the corresponding samples.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @importFrom jsonlite fromJSON
#' @importFrom SpatialFeatureExperiment addImg mirrorImg
#'
#' @export
load_visium_msfe <- function(dataset) {
  ## Update object name based on the dataset
  if (dataset == "lung") {
    object_name <- "msfe_2"
  } else if (dataset == "prostate") {
    object_name <- "msfe"
  } else {
    stop("Unsupported dataset. Please choose either 'lung' or 'prostate'.")
  }

  ## Use the `data` function to load the object
  data(object_name, package = "STExplorer")

  ## The object is now loaded into the environment, but we need to assign it to a variable
  msfe <- get(object_name, envir = .GlobalEnv)

  ## Get the sample names and directories
  sampleNames <- unlist(msfe@sample_ids)
  names(sampleNames) <- sampleNames  # Ensure sampleNames is named with its items

  ## Define the base directory for the sample data
  sampleDir_base <- system.file("extdata", "visium", package = "STExplorer")

  ## Define sampleDir based on the dataset
  if (dataset == "lung") {
    sampleDir <- list(
      Fibrotic = file.path(sampleDir_base, "lung", "Fibrotic"),
      Healthy  = file.path(sampleDir_base, "lung", "Healthy")
    )
  } else if (dataset == "prostate") {
    sampleDir <- list(H2_5 = file.path(sampleDir_base, "prostate", "H2_5"))
  }

  ## Process each sample
  for (id in sampleNames) {
    # Retrieve the SFE object for the current sample
    sfe <- getSFE(msfe, id)

    # Remove existing image data
    sfe@int_metadata$imgData <- NULL

    # Load the scale factors from the JSON file
    scaleF <- fromJSON(txt = file.path(sampleDir[[id]],
                                       "outs/spatial",
                                       "scalefactors_json.json"))

    # Add the low-resolution image to the SFE object
    sfe <- addImg(sfe,
                  file.path(sampleDir[[id]],
                            "outs/spatial/tissue_lowres_image.png"),
                  sample_id = id,
                  image_id = "lowres",
                  scale_fct = scaleF[["tissue_lowres_scalef"]])

    # Mirror the image for the SFE object
    sfe <- mirrorImg(sfe, sample_id = id, image_id = "lowres")

    # Add the modified SFE object back to the msfe object
    msfe <- addSFE(msfe, sfe, id)

    # Housekeeping: Remove the temporary SFE object from memory
    rm(sfe)
  }

  ## Assign the updated msfe object back to the global environment
  assign(object_name, msfe, envir = .GlobalEnv)
}

# ---------------------------------------------------------------------------- #
#  ########## INTERNAL FUNCTIONS ASSOCIATED WITH DISTANCE MATRIX ############
# ---------------------------------------------------------------------------- #
#' Internal Function: Add Distance Matrix to SpatialFeatureExperiment
#'
#' This function calculates the distance matrix for spatial coordinates and
#' adds it to the metadata of a SpatialFeatureExperiment object.
#'
#' @param sfe A SpatialFeatureExperiment object containing spatial coordinates.
#' @param p Numeric, the power parameter for the distance metric (default is 2
#' for Euclidean distance).
#' @param ... Additional arguments to be passed to the gw.dist function.
#'
#' @return A modified SpatialFeatureExperiment object with the distance matrix
#' added to its metadata.
#'
#' @details The function calculates the distance matrix for the spatial
#' coordinates in the SpatialFeatureExperiment object using the specified
#' distance metric (Euclidean, Manhattan, or Minkowski). The resulting distance
#' matrix is stored in the metadata of the SpatialFeatureExperiment object.
#'
#' @seealso \code{\link[GWmodel]{gw.dist}},
#' \code{\link[SpatialFeatureExperiment]{spatialCoords}},
#' \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords distance matrix, spatial coordinates, SpatialFeatureExperiment
#'
#' @rdname dot-int_addDistMat
#'
#' @importFrom GWmodel gw.dist
#'
.int_addDistMat <- function(sfe, p, ...) {
  ## Generate distance matrix
  dMat <- gw.dist(spatialCoords(sfe), p = p, ...)
  dimnames(dMat) <- list(NULL, colData(sfe)$Barcode)

  ## Determine the name of the distance metric
  if (p == 2) {
    dMetric <- "euclidean"
  } else if (p == 1) {
    dMetric <- "manhattan"
  } else if (p > 2) {
    dMetric <- paste0("minkowski_", p)
  }

  ## Store the distance matrix in the metadata of the SFE object
  sfe@metadata$dMat[[dMetric]] <- dMat

  return(sfe)
}

#' Internal Function: Get Distance Matrix from SpatialFeatureExperiment
#'
#' This function retrieves the distance matrix stored in the metadata of a
#' SpatialFeatureExperiment object.
#'
#' @param sfe A SpatialFeatureExperiment object containing a distance matrix in
#' its metadata.
#' @param dMetric The name of the distance metric used for the stored matrix
#' (Euclidean, Manhattan, or Minkowski).
#'
#' @return The distance matrix retrieved from the metadata of the
#' SpatialFeatureExperiment object.
#'
#' @details The function retrieves the distance matrix stored in the metadata
#' of the SpatialFeatureExperiment object. Users need to specify the distance
#' metric (Euclidean, Manhattan, or Minkowski) used for the stored matrix.
#'
#' @seealso \code{\link[SpatialFeatureExperiment]{colData}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords distance matrix, metadata, SpatialFeatureExperiment
#'
#' @rdname dot-int_getDistMat
#'
.int_getDistMat <- function(sfe, dMetric) {
  if (is.null(dMetric)) {
    dMetric <- names(sfe@metadata$dMat)[1]
  }

  dMat <- sfe@metadata$dMat[[dMetric]]

  return(dMat)
}



#' Internal: Generate a Gradient of Colours for a Single Colour
#'
#' This function creates a gradient of colours from a single colour to white,
#' with a specified number of steps.
#'
#' @param colour A character string representing the colour to use for the gradient.
#' @param steps An integer specifying the number of steps in the gradient.
#'
#' @return A vector of colour values representing the gradient.
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @rdname dot-int_getGradient
#'
#' @keywords internal
#'
.int_getGradient <- function(colour, steps) {
  # Convert colour to RGB values
  color_rgb <- grDevices::col2rgb(colour)

  # Calculate gradient steps
  steps2 <- steps + 2
  gradient_steps <- lapply(seq(0, 1, length.out = steps2), function(step) {
    grDevices::rgb(
      red = color_rgb[1] * step + 255 * (1 - step),
      green = color_rgb[2] * step + 255 * (1 - step),
      blue = color_rgb[3] * step + 255 * (1 - step),
      maxColorValue = 255
    )
  })

  # Return the gradient steps, excluding the first and last
  return(unlist(gradient_steps)[2:(steps + 1)])
}


# ---------------------------------------------------------------------------- #
#  ################# INTERNAL FUNCTIONS FOR MISC UTILITIES ###################
# ---------------------------------------------------------------------------- #
#' Internal: Format with leading zero
#'
#' This function appends the number `0` only for the numbers 1 to 9.
#'
#' @param x the number
#' @rdname dot-int_addLeadingZero
#' @author Eleftherios (Lefteris) Zormpas
#'
.int_addLeadingZero <- function(x) {
  if (x < 10) {
    return(paste0("0", x))
  } else {
    return(as.character(x))
  }
}


#' Internal: Find a file in a directory using regex
#'
#' This function returns the pathway for a file using regex
#'
#' @param directory character string. The directory.
#' @param pattern character string. The regex pattern.
#'
#' @rdname dot-int_findFileWithRegex
#' @author Eleftherios (Lefteris) Zormpas
#'
.int_findFileWithRegex <- function(directory, pattern) {
  # List all files in the directory
  files <- list.files(directory, full.names = TRUE)

  # Use grep to filter files that match the regex pattern
  matched_files <- grep(pattern, files, value = TRUE)

  # Check if any files matched
  if (length(matched_files) > 0) {
    # Normalize the path (optional, depending on your needs)
    matched_files <- normalizePath(matched_files)
    # Return the first matched file path
    return(matched_files[1])
  } else {
    # No file matched the pattern
    message("No file matching the pattern found.")
    return(NULL)
  }
}


#' INTERNAL: Get Sample IDs from a named List
#'
#' This internal function extracts sample IDs from a list, based on the
#' provided sample ID specification.
#'
#' @param list A list object containing sample information.
#' @param sample_id Either a logical vector indicating the samples to include
#' (if TRUE, all samples are included), or a character vector specifying the
#' sample IDs to include. It is suggested to use a character vector to specify
#' a specific sample, or NULL to select the first available sample.
#'
#' @return A character vector containing the selected sample IDs.
#'
#' @keywords internal
#'
#' @rdname dot-int_getListSmplIDs
.int_getListSmplIDs <- function(list, sample_id) {
  if (is.null(sample_id)) {
    ids <- names(list)[1]
  } else if (is.character(sample_id)) {
    ids <- sample_id
  } else if (sample_id) {
    ids <- names(list)
  }

  return(ids)
}


#' Internal Function: .int_getSmplIDs
#'
#' Description: Fetches sample IDs.
#'
#' @param sfe A SpatialFeatureExperiment object.
#' @param sample_id A character string, \code{TRUE}, or \code{NULL} specifying
#' sample/image identifier(s). If \code{TRUE}, all samples/images are
#' considered. If \code{NULL}, the first available entry is considered.
#'
#' @return Returns a character vector of sample IDs for plotting.
#'
#' @details This function fetches sample IDs based on the input parameters,
#' providing flexibility for customizing the sample IDs for plotting.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @keywords internal
#'
#' @rdname dot-int_getSmplIDs
#'
#' @importFrom DelayedArray unique
#'
.int_getSmplIDs <- function(sfe, sample_id = NULL) {
  ## Fetch the required sample IDs
  if (is.null(sample_id)) {
    sample_id <- DelayedArray::unique(colData(sfe)$sample_id)[1]
  } else if (is.character(sample_id)) {
    sample_id <- sample_id
  } else if (sample_id) {
    sample_id <- DelayedArray::unique(colData(sfe)$sample_id)
  }

  return(sample_id)
}
