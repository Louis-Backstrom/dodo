#' @title Strauss & Sadler's (1989) "Classical" Model
#'
#' @description
#' Equations 8 and 20 from Strauss & Sadler 1989. Estimates a one-tailed
#' \eqn{1 - \alpha} confidence interval and point estimate on the time of
#' extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#'
#' @returns a `list` object with the original parameters and the point estimate
#' and confidence interval included as elements. The confidence interval is a
#' two-element numeric vector called `conf.int`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Strauss, D., & Sadler, P. M. (1989). Classical Confidence Intervals and
#' Bayesian Probability Estimates for Ends of Local Taxon Ranges. Mathematical
#' Geology, 21(4), 411-421. \doi{10.1007/Bf00897326}
#'
#' **Other References**
#'
#' Marshall, C. R. (1990). Confidence intervals on stratigraphic ranges.
#' *Paleobiology*, 16(1), 1-10. \doi{10.1017/S0094837300009672}
#'
#' McFarlane, D. A. (1999). A Comparison of Methods for the Probabilistic
#' Determination of Vertebrate Extinction Chronologies. In R. D. E. MacPhee
#' (Ed.), *Extinctions in Near Time* (pp. 95-103). Springer US.
#' \doi{10.1007/978-1-4757-5202-1_5}
#'
#' Rivadeneira, M. M., Hunt, G., & Roy, K. (2009). The use of sighting records
#' to infer species extinctions: an evaluation of different methods. *Ecology*,
#' 90(5), 1291-1300. \doi{10.1890/08-0316.1}
#'
#' @export

SS89F1 <- function(records, alpha = 0.05) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # Calculate point estimate
  estimate <- (n * max(records) - min(records)) / (n - 1)

  # Calculate width of confidence interval
  x <- (alpha ^ (-1 / (n - 1)) - 1) * (max(records) - min(records))

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    estimate = estimate,
    conf.int = c(max(records), max(records) + x)
  )

  return(output)

}
