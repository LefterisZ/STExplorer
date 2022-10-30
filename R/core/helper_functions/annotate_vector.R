#' @name annotate_vector
#' 
#' @description a function to annotate any vector of of ENSGene IDs.
#' 
#' @param data A vector with an ENSGIDs as input.
#' @param biomart Biomart annotations of any organism. Can be generated using 
#'                the \code{read_biomart}.
#' @param add a vector of biomart columns to add to the input data frame. 
#'            Defaults to gene_name, biotype and description.
#' @param id position or name of column in biomart that matches the ENSGIDs 
#'           column or the rows names in input the data frame.
#' @param check.names argument for \code{data.frame} function at the end.
#' @param output takes one of two strings; "data.frame" (default) and "vector".
#'               The "data.frame" option will output the input vector alongside 
#'               the annotation columns in a data.frame. The "vector" option 
#'               will output a vector with a single annotation type for the 
#'               genes present in the input vector.
#' @param annot.col the annotation column to be returned as vector. Only needed 
#'                  if output = "vector"
#' 
#' @export

annotate_vector <- function(data, biomart, add, id = 1, check.names = FALSE,
                          output = "data.frame", annot.col = NULL){
    #### Match input vector with biomart rows ####
    if(!is.vector(data)){
        stop("Please check your input data. It must be a vector. If you wish to
             annotate a data.frame then use annotate_data.frame function.")
    } else {
        biom.dt <- match(data, biomart[[id]])
    }
    
    #### Check if user provided specific columns in add argument ####
    if(missing(add)){
        add <- c("gene_name", "biotype", "description") # If not provided, give default columns.
        if("human_homolog" %in% names(biomart)){ # If not human, biomart search for human homologs column
            add <- c(add, "human_homolog")
        }
    }
    
    #### Check if there was a successfull match between df and biomart ####
    if(all(is.na(biom.dt))){
        stop("ENSG IDs in data do not match column ", id, " in biomart table")
    }
    
    if(any(is.na(biom.dt))){
        message(sum(is.na(biom.dt)), " rows in data are missing from biomart table")
        warning("Rows in data missing from the biomart table, it might mean that
                you are using the wrong ensembl version to create the biomart.
                Please check the ensembl version you used for the annotation of 
                your data and generate a biomart table using the 'version' 
                parameter in the 'create_biomart' function")
        missing <- biomart[is.na(biom.dt)]
    }
    
    #### Output the data ####
    if(output == "data.frame"){
        out <- data.frame(id = data, # Combine to a new data frame
                          biomart[biom.dt,  add], 
                          stringsAsFactors = FALSE,
                          check.names = check.names)
    } else if(output == "vector") {
        out <- pull(biomart[biom.dt,  annot.col]) # Use pull from dplyr to get the values in a vector
    } else {
        stop("Please give a valid output value: 'data.frame' OR 'vector'")
    }
    
    return(out)

}
