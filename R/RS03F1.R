#' Roberts & Solow's (2003) "Optimal Linear Estimation" model
#'
#' @description
#' Equations from Roberts & Solow 2003 (Equations 12-18 from Solow 2005).
#' Estimates a p-value for testing competing hypotheses of
#' extinction/non-extinction, and a \eqn{1 - \alpha} confidence interval and
#' point estimate on the time of extinction.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param conf.int what kind of confidence interval to present. Valid options
#' are `"two-sided"` (the default) or `"one-sided"`.
#' @param k number of most recent sighting records to use (defaults to 10, or
#' the whole sighting record, if there are fewer than 10 records).
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
#' Clements, C. F., Collen, B., Blackburn, T. M., & Petchey, O. L. (2014).
#' Effects of recent environmental change on accuracy of inferences of
#' extinction status. *Conservation Biology*, 28(4), 971-981.
#' \doi{10.1111/cobi.12329}
#'
#' @examples
#' # Run the Dodo analysis from Roberts & Solow 2003
#' RS03F1(dodos, test.time = 2002)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' RS03F1(curlew$ccon, test.time = 2022)
#' }
#'
#' @export

RS03F1 <- function(records, alpha = 0.05, conf.int = "two-sided",
                   k = ifelse(length(records) >= 10, 10, length(records)),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {
  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # Set up function components
  sights <- sort(records, decreasing = TRUE)[1:k]
  v <- (1 / (k - 1)) *
    sum(log((sights[1] - sights[k]) / (sights[1] - sights[2:(k - 1)])))
  e <- matrix(rep(1, k), ncol = 1)
  if (conf.int == "two-sided") {
    SL <- (-log(1 - alpha / 2) / length(sights))^-v
    SU <- (-log(alpha / 2) / length(sights))^-v
  } else if (conf.int == "one-sided") {
    SL <- (-log(1 - 0) / length(sights))^-v
    SU <- (-log(alpha) / length(sights))^-v
  } else {
    stop("invalid value for conf.int, should be either one-sided or two-sided")
  }
  lambda <- outer(1:k, 1:k, myfun, v = v)
  lambda <- ifelse(lower.tri(lambda), lambda, t(lambda))
  a <- as.vector(solve(t(e) %*% solve(lambda) %*% e)) * solve(lambda) %*% e

  # Calculate model estimates
  conf.int.lower <- max(sights) + ((max(sights) - min(sights)) / (SL - 1))
  conf.int.upper <- max(sights) + ((max(sights) - min(sights)) / (SU - 1))
  estimate <- sum(t(a) %*% sights)
  p.value <- exp(-k * ((test.time - sights[1]) / (test.time - sights[k]))^
    (1 / v))

  if (conf.int.lower >= conf.int.upper) {
    warning("Confidence Interval estimation produced an invalid interval!")
  }

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

#' @title myfun from sExtinct package
#'
#' @description
#' Helper function. Modified myfun function from sExtinct package. For RS03F1().
#'
#' @param i a number.
#' @param j a number.
#' @param v a number.
#'
#' @returns a number.
#'
#' @references
#' **Key Reference**
#'
#' Clements, C. F., Collen, B., Blackburn, T. M., & Petchey, O. L. (2014).
#' Effects of recent environmental change on accuracy of inferences of
#' extinction status. *Conservation Biology*, 28(4), 971-981.
#' \doi{10.1111/cobi.12329}
#'
#' **Other References**
#'
#' Roberts, D. L., & Solow, A. R. (2003). Flightless birds: when did the dodo
#' become extinct? *Nature*, 426(6964), 245. \doi{10.1038/426245a}
#'
#' @noRd

myfun <- function(i, j, v) {
  return(exp((lgamma(2 * v + i) + lgamma(v + j)) - (lgamma(v + i) + lgamma(j))))
}
