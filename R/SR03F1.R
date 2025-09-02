#' Solow & Roberts' (2003) "Non-parametric" Model
#'
#' @description
#' Equation 4 from Solow & Roberts 2003 and Equations 9-11 from Solow 2005.
#' Estimates a p-value for testing competing hypotheses of
#' extinction/non-extinction, and a one-tailed \eqn{1 - \alpha} confidence
#' interval and point estimate on the time of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#'
#' @returns a `list` object with the original parameters and the p-value, point
#' estimate, and confidence interval included as elements. The confidence
#' interval is a two-element numeric vector called `conf.int`.
#'
#' @note
#' All sighting records are assumed to be certain.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Roberts, D. L. (2003). A nonparametric test for extinction
#' based on a sighting record. *Ecology*, 84(5), 1329-1332.
#' \doi{10.1890/0012-9658(2003)084[1329:ANTFEB]2.0.CO;2}
#'
#' **Other References**
#'
#' Robson, D. S., & Whitlock, J. H. (1964). Estimation of a truncation point.
#' *Biometrika*, 51(1-2), 33-39. \doi{10.1093/biomet/51.1-2.33}
#'
#' Solow, A. R. (2003). Estimation of stratigraphic ranges when fossil finds
#' are not randomly distributed. *Paleobiology*, 29(2), 181-185.
#' \doi{10.1666/0094-8373(2003)029<0181:EOSRWF>2.0.CO;2}
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#' *Mathematical Biosciences*, 195(1), 47-55. \doi{10.1016/j.mbs.2005.02.001}
#'
#' Rivadeneira, M. M., Hunt, G., & Roy, K. (2009). The use of sighting records
#' to infer species extinctions: an evaluation of different methods. *Ecology*,
#' 90(5), 1291-1300. \doi{10.1890/08-0316.1}
#'
#' @export

SR03F1 <- function(records, alpha = 0.05,
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # Calculate p-value
  p.value <- (records[n] - records[n - 1]) / (test.time - records[n - 1])

  # Calculate point estimate
  estimate <- records[n] + (records[n] - records[n - 1])

  # Calculate width of confidence interval
  x <- ((1 - alpha) / alpha) * (records[n] - records[n - 1])

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    test.time = test.time,
    p.value = p.value,
    estimate = estimate,
    conf.int = c(records[n], records[n] + x)
  )

  return(output)
}



