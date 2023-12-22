# Install and load the 'testthat' package
# install.packages("testthat")
library(testthat)

# Load the gearyw function and its optimized version
# source("gearyw.R")  # Replace with the actual file path

# Define a test case
test_that("calculate_gearyw produces the same result for both versions", {
  # Generate sample data
  set.seed(123)
  nb <- list(as.integer(c(2, 3)), as.integer(c(1, 3)), as.integer(c(1, 2)))
  weights <- list(c(0.5, 0.3), c(0.5, 0.2), c(0.3, 0.2))
  x <- c(1, 2, 3)
  card <- c(2, 2, 2)
  ftype <- TRUE
  zero.policy <- TRUE

  # Calculate results using the original and optimized versions
  result_original <- calculate_gearyw_original(nb, weights, x, card, ftype, zero.policy)
  # result_R <- calculate_gearyw_R(nb, weights, x, card, ftype, zero.policy)
  # result_R_optimized <- calculate_gearyw_optimized(nb, weights, x, card, ftype, zero.policy)
  result_R_optimized <- calculate_gearyw_optimized(nb, weights, x, card, ftype, zero.policy)

  # Check if the results are equal
  # expect_equal(result_original, result_R)
  #expect_equal(result_original, result_R_optimized)
  expect_equal(result_original, result_R_optimized)

})

# Define the original version of the function
calculate_gearyw_original <- function(listw_neighbours, listw_weights, x, card, ft, zero.policy) {
  .Call("gearyw", listw_neighbours, listw_weights,
        as.numeric(x), as.integer(card),
        as.logical(zero.policy), as.logical(ft), PACKAGE = "spdep")
}

# calculate_gearyw_R <- function(nb, weights, x, card, ftype, zero.policy) {
#   n <- length(card)
#   ans <- numeric(n)
#
#   for (i in seq_len(n)) {
#     if (card[i] == 0) {
#       if (zero.policy) {
#         ans[i] <- 0
#       } else {
#         ans[i] <- NA
#       }
#     } else {
#       sum_val <- 0
#       xi <- x[i]
#
#       for (j in seq_len(card[i])) {
#         k <- nb[[i]][j]
#         wt <- weights[[i]][j]
#         diff <- (xi - x[k])
#         if (ftype) {
#           res <- diff^2
#         } else {
#           res <- abs(diff)
#         }
#         sum_val <- sum_val + wt * res
#       }
#       ans[i] <- sum_val
#     }
#   }
#
#   return(ans)
# }

# Define the optimized version of the function
# calculate_gearyw_optimized <- function(nb, weights, x, card, ftype, zero.policy) {
#   n <- length(card)
#   ans <- numeric(n)
#
#   for (i in seq_len(n)) {
#     if (card[i] == 0) {
#       if (zero.policy) {
#         ans[i] <- 0
#       } else {
#         ans[i] <- NA
#       }
#     } else {
#       sum_val <- 0
#       xi <- x[i]
#
#       for (j in seq_len(card[i])) {
#         k <- nb[[i]][j]
#         wt <- weights[[i]][j]
#         diff <- (xi - x[k])
#         if (ftype) {
#           res <- diff^2
#         } else {
#           res <- abs(diff)
#         }
#         sum_val <- sum_val + wt * res
#       }
#       ans[i] <- sum_val
#     }
#   }
#
#   return(ans)
# }

# Define the optimized version of the function
calculate_gearyw_optimized <- function(nb, weights, x, card, ftype, zero.policy) {
  n <- length(card)
  ans <- numeric(n)

  for (i in seq_len(n)) {
    if (card[i] == 0) {
      if (zero.policy) {
        ans[i] <- 0
      } else {
        ans[i] <- NA
      }
    } else {
      xi <- x[i]
      diffs <- xi - nb[[i]]

      if (ftype) {
        res <- diffs^2
      } else {
        res <- abs(diffs)
      }

      sum_val <- sum(weights[[i]] * res)
      ans[i] <- sum_val
    }
  }

  return(ans)
}
