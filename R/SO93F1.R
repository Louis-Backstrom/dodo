#' @title Solow's (1993) "Classical" model
#'
#' @description
#' Equation 2 from Solow 1993 (Equation 3 from Solow 2005), and Equations 4
#' and 5 from Solow 2005. Estimates a p-value for testing competing hypotheses
#' of extinction/non-extinction, and a one-tailed \eqn{1 - \alpha} confidence
#' interval and point estimate on the time of extinction.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#'
#' @returns a `list` object with the original parameters and the p-value, point
#' estimate, and confidence interval included as elements. The confidence
#' interval is a two-element numeric vector called `conf.int`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (1993). Inferring Extinction from Sighting Data. *Ecology*,
#' 74(3), 962-964. \doi{10.2307/1940821}
#'
#' **Other References**
#'
#' Solow, A., & Helser, T. (2000). Detecting extinction in sighting data. In S.
#' Ferson & M. Burgman (Eds.), *Quantitative methods for conservation biology*
#' (pp. 1-6). Springer-Verlag. \doi{10.1007/b97704}
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#' *Mathematical Biosciences*, 195(1), 47-55. \doi{10.1016/j.mbs.2005.02.001}
#'
#' Rivadeneira, M. M., Hunt, G., & Roy, K. (2009). The use of sighting records
#' to infer species extinctions: an evaluation of different methods. *Ecology*,
#' 90(5), 1291-1300. \doi{10.1890/08-0316.1}
#'
#' @seealso [SO93B1()]
#'
#' @examples
#' # Run the Caribbean Monk Seal analysis from Solow 1993
#' SO93F1(monk_seal, test.time = 1992)
#' # Run an example analysis using the Slender-billed Curlew data
#' SO93F1(curlew$ccon, init.time = 1817, test.time = 2022)
#'
#' @export

SO93F1 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {
  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # If using first record as init.time, remove this from the record sequence
  if (init.time == min(records)) {
    records <- tail(records, -1)
    n <- n - 1
  }

  # Determine length of sighting record
  tn <- max(records) - init.time

  # Determine current time
  bigT <- test.time - init.time

  # Calculate p-value
  p.value <- (tn / bigT)^n

  # Calculate point estimate
  estimate <- ((n + 1) / n) * tn

  # Calculate relative width of confidence interval
  x <- alpha^(1 / n)

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.value = p.value,
    estimate = init.time + estimate,
    conf.int = c(max(records), init.time + tn / x)
  )

  return(output)
}
