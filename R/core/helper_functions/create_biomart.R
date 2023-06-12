#' @name create_biomart
#' 
#' @description a function that creates biomart tables containing ENSembl 
#'              annotations. It builds around the biomaRt package.
#' 
#' @param organism The first letter of genus and full species name like hsapiens.
#'                 A list of such names can be found in \code{listEnsembl}. A 
#'                 few common names are accepted for human, mouse, rat, 
#'                 zebrafish, fruitfly, pig, worm and yeast. 
#' @param attributes vector of column names to pass to \code{getBM}, default is:
#' ensembl_gene_id, external_gene_name, gene_biotype, chromosome_name,
#' start_position, end_position, strand, description and transcript_count
#' @param version Ensembl version for previous releases
#' @param patch Keep features on patches starting with CHR_, default FALSE
#' @param \dots additional options like filters and values passed to \code{getBM} or
#' \code{listAttributes}
#' 
#' @export

create_biomart <- function(organism = "human", attributes, version = NULL, 
                           patch = FALSE, ...){
    
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
        warning("Organism not in the common names list.\nPlease be sure you are using the right format:\nfirst letter of genus and full species name like hsapiens")
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
                                      version = version)
    
    #### Perform default search in the downloaded biomart ####
    if (missing(attributes)) {
        # Build default attributes list
        default.attr <- c("ensembl_gene_id","external_gene_name", "gene_biotype",
                          "chromosome_name", "start_position", "end_position",
                          "strand", "description", "transcript_count")
        # Get the biomart
        bm <- biomaRt::getBM(attributes = default.attr, mart = ensembl.bm, ...)
        
        # Replace long names like ensembl_gene_id with sorter ones
        names(bm)[1:6] <- c("id", "gene_name", "biotype", "chromosome", "start", "end")
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
                    message("Removed ", n.bm - nrow(bm), " features on patch CHR_*")
                } 
            }
        }
    }
    message("Downloaded ", nrow(bm), " features")
    
    
    # Return a data frame
    bm <- tibble::as_tibble(bm)
    bm
}
