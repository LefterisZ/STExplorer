#' Annotate a data frame
#'
#' @name annotateDataFrame
#'
#' @description
#' A function to annotate any data frame of data or biological data that has a
#' column of ENSGene IDs.
#'
#' @param data A data frame with an ENSGID column or rownames as input.
#' @param biomart Biomart annotations of any organism. Can be generated using
#' the \code{createBiomart} function.
#' @param add a vector of biomart columns to add to the input data frame.
#' Defaults to gene_name, biotype and description.
#' @param id position or name of column in biomart that matches the ENSGIDs
#' column or the rows names in input the data frame.
#' @param column IT MUST BE AN INTEGER: the number of the ENSGIDs column from
#' the input data frame. If left to NA then the function will look to match the
#' rownames. If no row names with ENSGIDs exist then returns an error.
#' @param check.names argument for \code{data.frame} function at the end.
#' @param output takes one of two strings; "data.frame" (default) and "vector".
#' The "data.frame" option will output the input data with the annotation
#' columns added. The "vector" option will output a vector with a single
#' annotation type for the genes present in the input data.
#' @param annot.col the annotation column to be returned as vector. Only needed
#' if output = "vector"
#'
#' @details
#' This function annotates a data frame (`data`) using a biomart database
#' (`biomart`). It matches the rows of the data frame with the corresponding
#' entries in the biomart based on the provided identifier (`id`) in the
#' biomart. If the data frame has a column with ENGene IDs then that specific
#' column (`column`) must be specified as an integer.
#'
#' If no specific biomart columns are provided in the `add` argument, the
#' function will add the following default columns: "gene_name", "biotype",
#' "description". If the biomart contains a column named "human_homolog", it
#' will be included in the annotation as well.
#'
#' If there is a successful match between the data frame and the biomart, the
#' annotated data is returned as a data frame or a vector based on the
#' specified `output` argument. The output data frame includes the original
#' data, the matching biomart columns, and an identifier column (`id`)
#' containing the row names of the input data frame.
#'
#' If some rows in the data frame are missing from the biomart table, a warning
#' message is displayed, indicating that the biomart might have been created
#' using a different Ensembl version. Users are advised to verify the Ensembl
#' version used for annotation and ensure consistency with the biomart table.
#'
#' @importFrom tibble as_tibble
#' @importFrom SpatialFeatureExperiment rowData
#'
#' @return An annotated data frame or vector based on the specified `output`
#' argument.
#'
#' @seealso \code{\link{createBiomart}}, \code{\link[tibble]{as_tibble}} to
#' convert the output to a tibble if necessary. \code{\link[biomaRt]{useMart}}
#' and \code{\link[biomaRt]{getBM}} to access and retrieve data from the biomart
#'  database.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' library(SpatialFeatureExperiment)
#'
#' # Load data
#' data(sfe)
#' df <- rowData(sfe)[,c(5:8)]
#'
#' # Create a biomart
#' biomart <- createBiomart(organism = "human")
#'
#' # Annotate tha dataframe
#' annotated_df <- annotateDataFrame(df, biomart = biomart)
#'
#' @export
annotateDataFrame <- function(data, biomart, add, id = 1, column = NA,
                              check.names = FALSE, output = "data.frame",
                              annot.col = NULL){
  ## Check if the input 'data' is a Data Frame
  if (!is.data.frame(data)) {
    stop("Error: Input 'data' is not a Data Frame.")
  }

  ## Match input df rows with biomart rows
  if (is.na(column)) {
    biom.dt <- match(rownames(data), biomart[[id]])
  } else if (is.numeric(column)) {
    biom.dt <- match(data[,column], biomart[[id]])
  } else {
    stop("Please check your data frame. No rownames with ENSGIDs were \n",
         "found to match the biomart database. If you provided an integer\n",
         "for a specific column please check again to ensure you provided\n",
         "the correct one.\n")
  }

  ## Check if user provided specific columns in add argument
  if (missing(add)) {
    add <- c("gene_name", "biotype", "description")
    if ("human_homolog" %in% names(biomart)) {
      # If not human, biomart search for human homologs column
      add <- c(add, "human_homolog")
    }
  }

  ## Check if there was a successfull match between df and biomart
  if (all(is.na(biom.dt)) && is.na(column)) {
    stop("Rownames in data do not match column ", id, " in biomart table.")
  } else if (all(is.na(biom.dt)) && is.integer(column)) {
    stop("ENSGIDs in provided column in data do not match column ", id,
         " in biomart table.")
  }

  if (any(is.na(biom.dt))) {
    message(sum(is.na(biom.dt)),
            " rows in data are missing from biomart table.")
    warning("Rows in data missing from the biomart table, it might ",
            "mean that you are using the wrong ensembl version to ",
            "create the biomart.\n",
            "Please check the ensembl version you used for the ",
            "annotation of your data when you mapped the reads to \n",
            "the genome and generate a biomart table using the ",
            "'version' parameter in the 'createBiomart' function.\n")
  }

  ## Output the data
  if (output == "data.frame") {
    out <- data.frame(id = rownames(data), # Combine to a new data frame
                      biomart[biom.dt,  add],
                      data,
                      stringsAsFactors = FALSE,
                      check.names = check.names)
  } else if (output == "vector") {
    out <- biomart[biom.dt,  annot.col]
  } else {
    stop("Please give a valid output value: 'data.frame' OR 'vector'")
  }

  tibble::as_tibble(out)
}


#' Annotate a vector
#'
#' @name annotateVector
#'
#' @description
#' A function to annotate any vector of ENSGene IDs.
#'
#' @param data A vector with ENSGIDs as input.
#' @param biomart Biomart annotations of any organism. Can be generated using
#' the \code{\link{createBiomart}} function from this package or the native
#' functions from the \code{biomaRt} package.
#' @param add A vector of biomart columns to add to the input data frame.
#' Defaults to "gene_name", "biotype", and "description".
#' @param id Position or name of the column in biomart that matches the ENSGIDs
#' column or the row names in the input data frame.
#' @param check.names Argument for \code{data.frame} function at the end.
#' @param output Takes one of two strings: "data.frame" (default) and "vector".
#' The "data.frame" option will output the input vector alongside the annotation
#' columns in a data.frame. The "vector" option will output a vector with a
#' single annotation type for the genes present in the input vector.
#' @param annot.col The annotation column to be returned as a vector. Only
#' needed if \code{output = "vector"}.
#'
#' @importFrom dplyr pull
#'
#' @seealso \code{\link{createBiomart}}
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' data <- c("ENSG00000139618", "ENSG00000141510", "ENSG00000198746")
#' biomart <- createBiomart(organism = "human")
#' annotated_vector <- annotateVector(data, biomart = biomart)
#'
#' @export
#'
annotateVector <- function(data, biomart, add, id = 1, check.names = FALSE,
                           output = "data.frame", annot.col = NULL){
  ## Match input vector with biomart rows
  if (!is.vector(data)) {
    stop("Please check your input data. It must be a vector. If you wish \n",
         "to annotate a data.frame then use `annotateDataFrame` function.\n")
  } else {
    biom.dt <- match(data, biomart[[id]])
  }

  ## Check if user provided specific columns in add argument
  if (missing(add)) {
    # If not provided, give default columns.
    add <- c("gene_name", "biotype", "description")
    if ("human_homolog" %in% names(biomart)) {
      # If not human, biomart search for human homologs column
      add <- c(add, "human_homolog")
    }
  }

  ## Check if there was a successful match between df and biomart
  if (all(is.na(biom.dt))) {
    stop("ENSG IDs in data do not match column ", id, " in biomart table")
  }

  if (any(is.na(biom.dt))) {
    message(sum(is.na(biom.dt)),
            " rows in data are missing from biomart table")
    warning("Rows in data missing from the biomart table, it might ",
            "mean that you are using the wrong ensembl version to ",
            "create the biomart.\n",
            "Please check the ensembl version you used for the ",
            "annotation of your data when you mapped the reads to \n",
            "the genome and generate a biomart table using the 'version' ",
            "parameter in the `createBiomart` function.\n")
  }

  ## Output the data
  if (output == "data.frame") {
    out <- data.frame(id = data, # Combine to a new data frame
                      biomart[biom.dt,  add],
                      stringsAsFactors = FALSE,
                      check.names = check.names)
  } else if (output == "vector") {
    if (is.null(annot.col)) {
      stop("You need to specify which column from the biomart you want as",
           "an output in vector format.\n",
           "  --> If you didn't supply your own columns with the 'add'\n",
           "      argument then the 'annot.col' argument needs to be one\n",
           "      of 'gene_name', 'biotype' or 'description'.\n",
           "  --> If you supplied your own columns with the 'add'\n",
           "      argument then make sure you typed the column of\n",
           "      interest correctly for the 'annot.col' argument.")
    }
    # Use pull from dplyr to get the values in a vector
    out <- pull(biomart[biom.dt,  annot.col])
  } else {
    stop("Please give a valid output value: 'data.frame' OR 'vector'")
  }

  return(out)

}


#' Create a BioMart table
#'
#' @name createBiomart
#'
#' @description A function that creates biomart tables containing Ensembl
#' annotations using the biomaRt package.
#'
#' @param organism The first letter of the genus and full species name, e.g.,
#' "hsapiens" for human. A list of common names can be found in
#' \code{listEnsembl}, or you can use accepted names for human, mouse, rat,
#' zebrafish, fruitfly, pig, worm, and yeast.
#' @param attributes A vector of column names to pass to \code{getBM}. The
#' default attributes include ensembl_gene_id, external_gene_name, gene_biotype,
#' chromosome_name, start_position, end_position, strand, description, and
#' transcript_count.
#' @param version Ensembl version for previous releases.
#' @param patch Keep features on patches starting with CHR_. Default is FALSE.
#' @param mirror Specify an Ensembl mirror to connect to. Valid options are
#' 'www', 'uswest', 'useast', 'asia'. If no mirror is specified, the primary
#' site at www.ensembl.org will be used. Mirrors are not available for the
#' Ensembl Genomes databases.
#' @param ... Additional options like filters and values passed to \code{getBM}
#' or \code{listAttributes}.
#'
#' @details The function downloads the biomart for the specified organism and
#' version, or the latest release if no version is provided. It performs a
#' default search using the attributes specified or the default attributes.
#' The resulting biomart table is returned as a data frame.
#'
#' @return A data frame containing the biomart table.
#'
#' @importFrom biomaRt listEnsembl useEnsembl getBM
#' @importFrom dplyr filter arrange
#' @importFrom tibble as_tibble
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @examples
#' # Create Human BioMart
#' bm <- createBiomart(organism = "human", version = NULL)
#' head(bm)
#'
#' # Create Mouse BioMart
#' bm <- createBiomart(organism = "mouse",
#' attributes = c("ensembl_gene_id","external_gene_name",
#' "gene_biotype", "chromosome_name"),
#' version = NULL)
#' head(bm)
#'
#' @seealso \code{\link{listEnsembl}}, \code{\link{useEnsembl}},
#' \code{\link{getBM}}
#'
#' @export
createBiomart <- function(organism = "human", attributes, version = NULL,
                          patch = FALSE, mirror = NULL, ...){

  ## Make a named vector with common organism names
  common.nms <- c(human     = "hsapiens",
                  mouse     = "mmusculus",
                  rat       = "rnorvegicus",
                  zebrafish = "drerio",
                  fruitfly  = "dmelanogaster",
                  pig       = "sscrofa",
                  worm      = "celegans",
                  yeast     = "scerevisiae")

  ## Check if organism is in the list
  if (organism %in% names(common.nms)) {
    organism <- common.nms[[organism]]
    message("BioMart organism used is: ", organism)
  } else if (organism %in% common.nms) {
    message("BioMart organism used is: ", organism)
  } else {
    warning(cat("Organism not in the common names list.\n",
                "Please be sure you are using the right format:\n",
                "first letter of genus and full species",
                "name like hsapiens\n"))
  }

  ## Prepare to search
  organism <- paste0(organism, "_gene_ensembl")
  release <- version
  if (is.null(version)) {
    x.bm <- biomaRt::listEnsembl()
    release <- x.bm$version[x.bm$biomart == "genes"]
    release <- gsub("Ensembl Genes ", "", release)
    message("Using latest Ensembl release ", release)
  } else {
    x.bm <- biomaRt::listEnsembl(version = version)
    message("Using requested Ensembl release ", release)
  }

  ## Get the organism specific ensembl biomart
  ensembl.bm <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = organism,
                                    version = version,
                                    mirror = mirror)

  ## Perform default search in the downloaded biomart
  if (missing(attributes)) {
    # Build default attributes list
    default.attr <- c("ensembl_gene_id","external_gene_name",
                      "gene_biotype", "chromosome_name", "start_position",
                      "end_position", "strand", "description",
                      "transcript_count")
    # Get the biomart
    bm <- biomaRt::getBM(attributes = default.attr, mart = ensembl.bm, ...)

    # Replace long names like ensembl_gene_id with sorter ones
    names(bm)[1:6] <- c("id", "gene_name", "biotype",
                        "chromosome", "start", "end")
    n <- length(unique(bm$id))

    # Drop source from description: [Source:MGI Symbol;Acc:MGI:102478]
    bm$description <- gsub(" \\[.*\\]$", "" , bm$description)

    # Remove white space in version 92
    if (release == 92) {
      bm$description <- trimws(bm$description)
      bm <- dplyr::arrange(bm, .data$id)
    }

    # Remove features starting with CHR_
    if (!patch) {
      n.bm <- nrow(bm)
      bm <- dplyr::filter(bm, substr(.data$chromosome,1,4) != "CHR_")
      if (n.bm != nrow(bm)) {
        message("Removed ", n.bm - nrow(bm),
                " features on patch CHR_*")
      }
    }

  } else {
    # Get the biomart
    bm <- biomaRt::getBM(attributes = attributes, mart = ensembl.bm, ...)

    # Remove features starting with CHR_
    if (!patch) {
      if ("chromosome_name" %in% colnames(bm)) {
        n.bm <- nrow(bm)
        bm <- dplyr::filter(bm,
                            substr(.data$chromosome_name,1,4) != "CHR_")
        if (n.bm != nrow(bm)) {
          message(cat("Removed ", n.bm - nrow(bm),
                      " features on patch CHR_*\n"))
        }
      }
    }
  }
  message("Downloaded ", nrow(bm), " features")


  # Return a data frame
  bm <- tibble::as_tibble(bm)
  bm
}


# ---------------------------------------------------------------------------- #
# # INTERNAL FUNCTIONS ASSOCIATED ANNOTATING/MATCHING GENENAMES <-> ENSGIDs ----
# ---------------------------------------------------------------------------- #
#' Internal Function: match gene names to EnsgIDs
#'
#' This internal function takes an SFE object and a vector of gene names or
#' EnsgIDs, checks if there are gene names and returns their EnsgIDs from the
#' SFE object.
#'
#' @param sfe An SFE object.
#' @param genes A vector of EnsgIDs or gene names. A mixture of both throws an
#' error.
#'
#' @returns a vector of EnsgIDs
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_matchNamesToEnsgID
#'
.int_matchNamesToEnsgID <- function(sfe, genes) {
  ## Validate gene input
  isEnsgID <- grepl("ENS", genes)
  if (sum(isEnsgID) == length(genes)) {
    inputType <- "all_ensg"
  } else if (sum(isEnsgID) == 0) {
    inputType <- "no_ensg"
    # assuming all are gene names, check that are present
    isExisting <- genes %in% rowData(sfe)[, "gene_name"]
    allExisting <- isExisting == length(genes)
  } else {
    inputType <- "mixed"
  }

  ## Handle potential errors
  ### A. Gene names not present in the SFE
  if (inputType == "no_ensg") {
    if (!allExisting) {
      stop("Gene names not found in the SFE object: ",
           paste(genes[!isExisting], collapse = " "))
    }
  }

  ### B. Mix of EnsgIDs and Gene Names in the input
  if (inputType == "mixed") {
    stop("Check your genes input!\n",
         "--> You have probably provided a mixture of ENSG IDs and gene names.",
         "\n Please provide either ENSG IDs OR gene names. NOT both.")
  }

  ## Match gene names based on the validated input
  if (inputType == "all_ensg") {
    return(genes)
  } else if (inputType == "no_ensg") {
    match <- rowData(sfe)[, "gene_name"] %in% genes
    return(rownames(sfe)[match])
  }
}

#' Internal Function: match EnsgIDs to gene names without the use of a biomart
#'
#' This function uses the SFE object to translate a vector of ENSGIDs into gene
#' names.
#'
#' @inheritParams .int_matchNamesToEnsgID
#' @param genes A vector of EnsgIDs.
#'
#' @returns a vector of Gene Names
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
#' @rdname dot-int_EnsgIDtoName
#'
.int_EnsgIDtoName <- function(sfe, genes) {
  ## Check columns exist and if not, return a meaningful error
  isExisting <- c("id", "gene_name") %in% colnames(rowData(sfe))
  if (sum(isExisting) < 2) {
    stop("The provided SFE object is missing column/s ",
         c("id", "gene_name")[isExisting],
         " from rowData.")
  }
  ## Find the EnsgID positions in the rowData and export alongside gene names
  selection <- rowData(sfe)[["id"]] %in% genes
  df <- rowData(sfe)[selection, c("id", "gene_name")] %>%
    as.data.frame()

  ## Left join the export with the genes vector
  out <- genes %>%
    as.data.frame() %>%
    rename("id" = ".") %>%
    left_join(df, by = "id") %>%
    select("gene_name")

  ## Return
  ## It looks weird but tried as.vector and it returns a list maybe because out
  ## is a list? See as.vector help page.
  return(out[["gene_name"]])
}
