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
    #### Match input vector with biomart rows ####
    if (!is.vector(data)) {
      cat("Please check your input data. It must be a vector. If you wish to\n",
          "annotate a data.frame then use annotate_data.frame function.\n")
      stop()
    } else {
        biom.dt <- match(data, biomart[[id]])
    }

    #### Check if user provided specific columns in add argument ####
    if (missing(add)) {
      # If not provided, give default columns.
        add <- c("gene_name", "biotype", "description")
        if ("human_homolog" %in% names(biomart)) {
          # If not human, biomart search for human homologs column
            add <- c(add, "human_homolog")
        }
    }

    #### Check if there was a successfull match between df and biomart ####
    if (all(is.na(biom.dt))) {
        stop("ENSG IDs in data do not match column ", id, " in biomart table")
    }

    if (any(is.na(biom.dt))) {
      message(sum(is.na(biom.dt)),
              " rows in data are missing from biomart table")
      warning(cat("Rows in data missing from the biomart table, it might",
                  "mean that you are using the wrong ensembl version to",
                  "create the biomart.\n",
                  "Please check the ensembl version you used for the",
                  "annotation of your data when you mapped the reads to the",
                  "genome and generate a biomart table using the 'version'",
                  "parameter in the 'create_biomart' function.\n"))
        missing <- biomart[is.na(biom.dt)]
    }

    #### Output the data ####
    if (output == "data.frame") {
        out <- data.frame(id = data, # Combine to a new data frame
                          biomart[biom.dt,  add],
                          stringsAsFactors = FALSE,
                          check.names = check.names)
    } else if (output == "vector") {
        if (is.null(annot.col)) {
          cat("You need to specify which column from the biomart you want as",
              "an output in vector format.\n",
              "  --> If you didn't supply your own columns with the 'add'\n",
              "      argument then the 'annot.col' argument needs to be one\n",
              "      of 'gene_name', 'biotype' or 'description'.\n",
              "  --> If you supplied your own columns with the 'add'\n",
              "      argument then make sure you typed the column of\n",
              "      interest correctly for the 'annot.col' argument.")
          stop()
        }
      # Use pull from dplyr to get the values in a vector
        out <- pull(biomart[biom.dt,  annot.col])
    } else {
        stop("Please give a valid output value: 'data.frame' OR 'vector'")
    }

    return(out)

}
