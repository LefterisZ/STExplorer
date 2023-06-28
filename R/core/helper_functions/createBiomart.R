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

    #### Make a named vector with common organism names ####
    common.nms <- c(human     = "hsapiens",
                    mouse     = "mmusculus",
                    rat       = "rnorvegicus",
                    zebrafish = "drerio",
                    fruitfly  = "dmelanogaster",
                    pig       = "sscrofa",
                    worm      = "celegans",
                    yeast     = "scerevisiae")

    #### Check if organism is in the list ####
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

    #### Prepare to search ####
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

    #### Get the organism specific ensembl biomart ####
    ensembl.bm <- biomaRt::useEnsembl(biomart = "ensembl",
                                      dataset = organism,
                                      version = version,
                                      mirror = mirror)

    #### Perform default search in the downloaded biomart ####
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
            bm <- dplyr::arrange(bm, id)
        }

        # Remove features starting with CHR_
        if (!patch) {
            n.bm <- nrow(bm)
            bm <- dplyr::filter(bm, substr(chromosome,1,4) != "CHR_")
            if (n.bm != nrow(bm)) {
                message("Removed ", n.bm - nrow(bm), " features on patch CHR_*")
            }
        }

    } else {
        # Get the biomart
        bm <- biomaRt::getBM(attributes = attributes, mart = ensembl.bm, ...)

        # Remove features starting with CHR_
        if (!patch) {
            if ("chromosome_name" %in% colnames(bm)) {
                n.bm <- nrow(bm)
                bm <- dplyr::filter(bm, substr(chromosome_name,1,4) != "CHR_")
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
