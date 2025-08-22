#' @title F(x) from Solow 2005
#'
#' @description
#' Helper function. Modified Equation 7 from Solow 2005.
#'
#' @param x generally \eqn{t_n} or \eqn{T}.
#' @param s sum of relative sighting times.
#' @param n number of sighting records.
#'
#' @returns a number.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#' *Mathematical Biosciences*, 195(1), 47-55. \doi{10.1016/j.mbs.2005.02.001}
#'
#' @noRd

Fx <- function(x, s, n) {
  is <- 1:floor(s / x)

  part1 <- (-1) ^ (is - 1)
  part2 <- choose(n, is)
  part3 <- (1 - (is * x / s)) ^ (n - 1)

  result <- 1 - sum(part1 * part2 * part3)
  rm(is, part1, part2, part3)

  return(result)

}

#' @title Falling Factorial
#'
#' @description
#' Calculates the falling factorial function.
#'
#' @param x an integer.
#' @param j an integer.
#'
#' @returns a number.
#'
#' @noRd

ffactorial <- function(x, j) {
  return(factorial(x) / factorial(x - j))
}
