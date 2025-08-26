#' @title McFarlane's (1999) "Median Gap" Model
#'
#' @description
#' Equation 3 from McFarlane 1999. Estimates a one-tailed \eqn{1 - \alpha}
#' confidence interval on the time of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#'
#' @returns a `list` object with the original parameters and the confidence
#' interval included as elements. The confidence interval is a two-element
#' numeric vector called `conf.int`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' McFarlane, D. A. (1999). A Comparison of Methods for the Probabilistic
#' Determination of Vertebrate Extinction Chronologies. In R. D. E. MacPhee
#' (Ed.), *Extinctions in Near Time* (pp. 95-103). Springer US.
#' \doi{10.1007/978-1-4757-5202-1_5}
#'
#' @export

MC99F1 <- function(records, alpha = 0.05) {

  # Sort records
  records <- sort(records)

  # Calculate median gap length
  i <- median(diff(sort(records)))

  # Calculate critical value
  crit <- log(alpha, base = 0.5)

  # Calculate width of confidence interval
  x <- i * crit

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    conf.int = c(max(records), max(records) + x)
  )

  return(output)

}
