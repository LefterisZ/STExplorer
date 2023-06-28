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
    #### Match input df rows with biomart rows ####
    if (is.na(column)) {
        biom.dt <- match(rownames(data), biomart[[id]])
    } else if (is.numeric(column)) {
        biom.dt <- match(data[,column], biomart[[id]])
    } else {
        cat("Please check your data frame. No rownames with ENSGIDs were \n",
            "found to match the biomart database. If you provided an integer\n",
            "for a specific column please check again to ensure you provided\n",
            "the correct one.\n")
        stop()
    }

    #### Check if user provided specific columns in add argument ####
    if (missing(add)) {
        add <- c("gene_name", "biotype", "description")
        if ("human_homolog" %in% names(biomart)) {
          # If not human, biomart search for human homologs column
            add <- c(add, "human_homolog")
        }
    }

    #### Check if there was a successfull match between df and biomart ####
    if (all(is.na(biom.dt)) && is.na(column)) {
        stop("Rownames in data do not match column ", id, " in biomart table")
    } else if (all(is.na(biom.dt)) && is.integer(column)) {
        cat("ENSGIDs in provided column in data do not match column ", id,
            " in biomart table")
        stop()
    }

    if (any(is.na(biom.dt))) {
        message(sum(is.na(biom.dt)),
                " rows in data are missing from biomart table")
        warning(cat("Rows in data missing from the biomart table, it might",
                    "mean that you are using the wrong ensembl version to",
                    "create the biomart.\n",
                    "Please check the ensembl version you used for the",
                    "annotation of your data when you mapped the reads to \n",
                    "the genome and generate a biomart table using the ",
                    "'version' parameter in the 'createBiomart' function.\n"))
    }

    #### Output the data ####
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
