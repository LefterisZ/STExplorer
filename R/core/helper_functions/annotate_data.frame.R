#' @name annotate_data.frame
#' 
#' @description a function to annotate any data frame of data or biological 
#'              data that has a column of ENSGene IDs.
#' 
#' @param data A data frame with an ENSGID column or rownames as input.
#' @param biomart Biomart annotations of any organism. Can be generated using 
#'                the \code{read_biomart}.
#' @param add a vector of biomart columns to add to the input data drame. 
#'            Defaults to gene_name, biotype and description.
#' @param id position or name of column in biomart that matches the ENSGIDs 
#'           column or the rows names in input the data frame.
#' @param column IT MUST BE AN INTEGER: the number of the ENSGIDs column from the 
#'               input data frame. If left to NA then the function will look to
#'               match the rownames. If no row names with ENSGIDs exist then 
#'               returns an error.
#' @param check.names argument for \code{data.frame} function at the end.
#' @param output takes one of two strings; "data.frame" (default) and "vector".
#'               The "data.frame" option will output the input data with the 
#'               annotation columns added. The "vector" option will output a 
#'               vector with a single annotation type for the genes present
#'               in the input data.
#' @param annot.col the annotation column to be returned as vector. Only needed 
#'                  if output = "vector"
#' 
#' @export

annotate_data.frame <- function(data, biomart, add, id = 1, column = NA, 
                                check.names = FALSE, output = "data.frame", 
                                annot.col = NULL){
    #### Match input df rows with biomart rows ####
    if(is.na(column)){
        biom.dt <- match(rownames(data), biomart[[id]])
    } else if(is.numeric(column)){
        biom.dt <- match(data[,column], biomart[[id]])
    } else {
        stop("Please check your data frame. No rownames with ENSGIDs were found
             to match the biomart database. If you provided an integer for a 
             specific column please check again to ensure you provided the 
             correct one.")
    }
    
    #### Check if user provided specific columns in add argument ####
    if(missing(add)){
        add <- c("gene_name", "biotype", "description")
        if("human_homolog" %in% names(biomart)){ # If not human, biomart search for human homologs column
            add <- c(add, "human_homolog")
        }
    }
    
    #### Check if there was a successfull match between df and biomart ####
    if(all(is.na(biom.dt)) && is.na(column)){
        stop("Rownames in data do not match column ", id, " in biomart table")
    } else if(all(is.na(biom.dt)) && is.integer(column)){
        stop("ENSGIDs in provided column in data do not match column ", id, " in biomart table")
    }
    
    if(any(is.na(biom.dt))){
        message(sum(is.na(biom.dt)), " rows in data are missing from biomart table")
        warning("Rows in data missing from the biomart table, it might mean that
                you are using the wrong ensembl version to create the biomart.
                Please check the ensembl version you used for the annotation of 
                your data and generate a biomart table using the 'version' 
                parameter in the 'create_biomart' function")
    }
    
    #### Output the data ####
    if(output == "data.frame"){
        out <- data.frame(id = rownames(data), # Combine to a new data frame
                          biomart[biom.dt,  add], 
                          data, 
                          stringsAsFactors = FALSE,
                          check.names = check.names)
    } else if(output == "vector") {
        out <- biomart[biom.dt,  annot.col]
    } else {
        stop("Please give a valid output value: 'data.frame' OR 'vector'")
    }
    
    tibble::as_tibble(out)
}
