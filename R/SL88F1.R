#' @title Springer & Lilje's (1988) "Broken Stick" model
#'
#' @description
#' Equation 3 from Springer & Lilje 1988. Estimates a one-tailed
#' \eqn{1 - \alpha} confidence interval on the time of extinction.
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
#' Springer, M., & Lilje, A. (1988). Biostratigraphy and Gap Analysis - the
#' Expected Sequence of Biostratigraphic Events. *Journal of Geology*, 96(2),
#' 228-236. \doi{10.1086/629212}
#'
#'**Other References**
#'
#' McFarlane, D. A. (1999). A Comparison of Methods for the Probabilistic
#' Determination of Vertebrate Extinction Chronologies. In R. D. E. MacPhee
#' (Ed.), *Extinctions in Near Time* (pp. 95-103). Springer US.
#' \doi{10.1007/978-1-4757-5202-1_5}
#'
#' Strauss, D., & Sadler, P. M. (1989). Classical Confidence Intervals and
#' Bayesian Probability Estimates for Ends of Local Taxon Ranges.
#' *Mathematical Geology*, 21(4), 411-421. \doi{10.1007/Bf00897326}
#'
#' @export

SL88F1 <- function(records, alpha = 0.05) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # Calculate sighting rate
  m <- (n - 1) / (max(records) - min(records))

  # Calculate width of confidence interval
  x <- log(alpha) / (-m)

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    conf.int = c(max(records), max(records) + x)
  )

  return(output)

}
