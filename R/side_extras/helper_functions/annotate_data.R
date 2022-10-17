#' @name annotate_data
#' 
#' @description a function to annotate any data frame of datas or biological 
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
#' @param check.names argument from \code{data.frame} function at the end.
#' 
#' @export

annotate_data <- function(data, biomart, add, id = 1, column = NA, check.names){
    #### Match input df rows with biomart rows ####
    if(is.na(column)){
        biom.dt <- match(rownames(data), biomart[[id]])
    } else if(is.integer(column)){
        biom.dt <- match(data[,column], biomart[[id]])
    } else {
        stop("Please check your data frame. No rownames with ENSGIDs were found.\n
             to match the biomart database. If you provided an integer for a \n
             specific column please check again to ensure you provided the \n
             correct one.")
    }
    
    #### Check if user provided specific columns in add argument ####
    if(missing(add)){
        add <- c("gene_name", "biotype", "description")
        if("human_homolog" %in% names(biomart)){ # If not human biomart search for human homologs column
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
    }
    
    #### Combine to a new data frame ####
    out <- data.frame(id = rownames(data), 
                      biomart[biom.dt,  add], 
                      data, 
                      stringsAsFactors = FALSE,
                      check.names = check.names)
    tibble::as_tibble(out)
}
