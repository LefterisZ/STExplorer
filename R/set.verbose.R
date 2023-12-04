#' Set up verbose handlers
#'
#' @name set.verbose
#'
#' @description
#' IMPORTANT: this function is not exported to be used alone.
#' It is used internally to set the requirements for verbosity using the
#' progressr package.
#'
#' @keywords internal
#'
#' @param .global Logical value indicating whether to set the global handler for
#' verbosity. If TRUE, the global handler will be set; if FALSE, it will not
#' be set.
#' @param .handlers A character vector specifying the types of handlers to be
#' set for verbosity. Possible values include "progress" and "beepr".
#'
#' @importFrom progressr handlers
#'
#' @details
#' This function sets up the requirements for verbosity through the progressr
#' package. By default, it sets the global handler for verbosity, which affects
#' all verbosity calls. Additionally, it allows specifying specific types of
#' handlers to be used for verbosity, such as "progress" and "beepr". The
#' function helps control the verbosity settings in the code.
#'
#' @author Eleftherios (Lefteris) Zormpas
#'
set.verbose <- function(.global = TRUE,
                        .handlers = c("progress", "beepr")){
  ## Set global handler?
  progressr::handlers(global = .global)
  ## Set progress types
  progressr::handlers(.handlers)
}
