#' @title Marshall's (1994) "Distribution-free" model
#'
#' @description
#' Equations 1 and 2 from Marshall 1994. Estimates a one-sided
#' \eqn{1 - \alpha} confidence interval on time of extinction; the upper bound
#' of this confidence interval is in itself a \eqn{1 - 2\gamma} confidence
#' interval with a lower and upper bound.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param gamma desired confidence probability of the bounds of the confidence
#' interval.
#'
#' @returns a `list` object with the original parameters and the p-value, point
#' estimate, and confidence interval included as elements. The confidence
#' interval is a three-element numeric vector called `conf.int`; the first
#' element is the lower bound of the confidence interval (the time of the last
#' record), the second element is the lower bound of the upper bound of the
#' confidence interval, and the third element is the upper bound of the upper
#' bound of the confidence interval.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant. The upper bound of the upper confidence interval is likely to
#' be NA, unless there are sufficient sightings or \eqn{\alpha} and \eqn{\gamma}
#' are sufficiently permissive.
#'
#' @references
#' **Key Reference**
#'
#' Marshall, C. R. (1994). Confidence intervals on stratigraphic ranges:
#' partial relaxation of the assumption of randomly distributed fossil
#' horizons. *Paleobiology*, 20(4), 459-469. \doi{10.1017/S0094837300012938}
#'
#' @examples
#' # Run the *Metrarabdotos* n. sp. 5 analysis from Marshall 1994
#' MA94F1(metrarabdotos, alpha = 0.5, gamma = 0.025)
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' MA94F1(curlew$ccon, alpha = 0.05, gamma = 0.05)
#' }
#'
#' @export

MA94F1 <- function(records, alpha = 0.05, gamma = 0.1) {
  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)
  N <- n - 1

  i <- 0:N

  test <- c()

  for (xi in i) {
    x <- 0:xi
    test[xi + 1] <- sum(choose(N, x) * (1 - alpha)^x * alpha^(N - x))
  }

  x_lower <- i[max(which(TRUE == (gamma > test)))]
  x_upper <- i[min(which(TRUE == ((1 - gamma) < test)))]
  if (x_upper == N) {
    x_upper <- NA
  }

  output <- list(
    records = records,
    alpha = alpha,
    gamma = gamma,
    conf.int = c(
      max(records),
      max(records) + sort(diff(records))[x_lower],
      max(records) + unique(sort(diff(records))[x_upper])
    ) # unique() to fix NAs
  )

  return(output)
}
