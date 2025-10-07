#' @title McFarlane's (1999) "Median Gap" model
#'
#' @description
#' Equation 3 from McFarlane 1999. Estimates a p-value for testing competing
#' hypotheses of extinction/non-extinction,and a one-tailed \eqn{1 - \alpha}
#' confidence interval on the time of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#'
#' @returns a `list` object with the original parameters and the p-value and
#' confidence interval included as elements. The confidence interval is
#' a two-element numeric vector called `conf.int`.
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
#' @examples
#' # Run an example analysis using the Woolly Mammoth data
#' MC99F1(mammoth, test.time = -11000)
#'
#' @export

MC99F1 <- function(records, alpha = 0.05,
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {
  # Sort records
  records <- sort(records)

  # Calculate median gap length
  i <- median(diff(records))

  # Calculate critical value
  crit <- log(alpha, base = 0.5)

  # Calculate width of confidence interval
  x <- i * crit

  # Calculate width of test gap
  w <- test.time - max(records)

  # Calculate p-value
  p.value <- exp((w / i) * log(0.5))

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    p.value = p.value,
    conf.int = c(max(records), max(records) + x)
  )

  return(output)
}
