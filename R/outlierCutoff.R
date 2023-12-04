#' Get outlier cut off value
#'
#' @name outlierCutoff
#' @description
#' This function identifies outliers in a vector of data based on a coefficient
#' and the interquartile range (IQR).
#'
#' @param dt A numeric vector of data.
#' @param coef The coefficient used to determine the cutoff for outliers.
#' The default value is 1.5, which is commonly used in boxplots to define
#' outliers as values that are more than 1.5 times the IQR away from the upper
#' or lower quartiles.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{out_up}: The upper cutoff value for outliers.
#'   \item \code{out_down}: The lower cutoff value for outliers.
#'   \item \code{outliers}: A logical vector indicating whether each data point
#'   is an outlier or not. \code{TRUE} indicates an outlier, and \code{FALSE}
#'   indicates a non-outlier.
#' }
#'
#' @importFrom stats quantile
#'
#' @examples
#' # Generate a vector of random data
#' data <- rnorm(100)
#'
#' # Detect outliers using the default coefficient (1.5)
#' outlierCutoff(data)
#'
#' # Detect outliers using a custom coefficient (2.0)
#' outlierCutoff(data, coef = 2.0)
#'
#' @export
outlierCutoff <- function(dt, coef = 1.5) {
  ## Calculate the quartiles of the data
  quantiles <- quantile(dt, probs = c(0.25, 0.75))

  ## Calculate the interquartile range (IQR)
  IQR <- quantiles[2] - quantiles[1]

  ## Calculate the cutoff values for outliers
  out_up <- quantiles[2] + coef * IQR
  out_down <- quantiles[1] - coef * IQR

  ## Identify outliers based on the cutoff values
  outliers <- dt < out_down | dt > out_up

  ## Create a list with the results
  res <- list(out_up = out_up, out_down = out_down, outliers = outliers)

  return(res)
}
