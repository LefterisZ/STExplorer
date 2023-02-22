#' @name set.verbose
#' 
#' @description a function to set the requirements for verbose through progressr
#'              package.
#' @param .global progressr package arguments for progressr::handlers. Used only
#'                when verbose = TRUE.
#' @param .handlers progressr package arguments for progressr::handlers             
#'


set.verbose <- function(.global = TRUE, 
                        .handlers = c("progress", "beepr")){
  ## Set global handler?
  handlers(global = .global)
  ## Set progress types
  handlers(.handlers)
}



