#' Roberts & Solow's (2003) "Optimal Linear Estimation" model
#'
#' @description
#' Equations from Roberts & Solow 2003 (Equations 12-18 from Solow 2005).
#' Estimates a p-value for testing competing hypotheses of
#' extinction/non-extinction, and a two-tailed \eqn{1 - \alpha} confidence
#' interval and point estimate on the time of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param k number of most recent sighting records to use (defaults to 5).
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
#' Roberts, D. L., & Solow, A. R. (2003). Flightless birds: when did the dodo
#' become extinct? *Nature*, 426(6964), 245. \doi{10.1038/426245a}
#'
#' **Other References**
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#' *Mathematical Biosciences*, 195(1), 47-55. \doi{10.1016/j.mbs.2005.02.001}
#'
#' Rivadeneira, M. M., Hunt, G., & Roy, K. (2009). The use of sighting records
#' to infer species extinctions: an evaluation of different methods. *Ecology*,
#' 90(5), 1291-1300. \doi{10.1890/08-0316.1}
#'
#' @export

RS03F1 <- function(records, alpha = 0.05, k = 10,
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # Set up function components
  e <- rep(1, k)
  nuhat <- (1 / (k - 1)) * sum(log(
    (records[n] - records[n - k + 1]) /
      (records[n] - records[(1:(k - 2)) + 1])))
  Lambda <- matrix(ncol = k, nrow = k)
  for (j in 1:k) {
    for (i in 1:k) {
      if (j <= i) {
        Lambda[i, j] <- (gamma(2 * nuhat + i) * gamma(nuhat + j)) /
          (gamma(nuhat + i) * gamma(j))
      } else {
        Lambda[i, j] <- Lambda[j, i]
      }
    }
  }
  w <- as.vector(solve(t(e) %*% solve(Lambda) %*% matrix(e))) *
    solve(Lambda) %*% matrix(e)
  c <- (-k / log(alpha)) ^ (-nuhat)

  # Calculate p-value
  p.value <- exp(-k * ((test.time - records[n]) /
                   (test.time - records[n - k + 1])) ^ (1 / nuhat))

  # Calculate point estimate
  estimate <- sum(w * records[n - 1:k + 1])

  # Calculate lower bound of confidence interval
  SL <- (-log(1 - alpha / 2) / k) ^ (-nuhat)
  conf.int.lower <- records[n] + (records[n] - records[n - k + 1]) / (SL - 1)

  # Calculate upper bound of confidence interval
  SU <- (-log(alpha / 2) / k) ^ (-nuhat)
  conf.int.upper <- records[n] + (records[n] - records[n - k + 1]) / (SU - 1)

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    k = k,
    test.time = test.time,
    p.value = p.value,
    estimate = estimate,
    conf.int = c(conf.int.lower, conf.int.upper)
  )

  return(output)

}


