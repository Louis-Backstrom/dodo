#' @title Solow's (1993) "Declining" Model
#'
#' @description
#' Equations 3 and 4 from Solow 1993 (Equation 7 from Solow 2005), and
#' Equation 8 from Solow 2005. Estimates a p-value for testing competing
#' hypotheses of extinction/non-extinction, and a one-tailed \eqn{1 - \alpha}
#' confidence interval and point estimate on the time of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case  this sighting is removed from the record.
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
#' Solow, A. R. (1993). Inferring Extinction in a Declining Population.
#' *Journal of Mathematical Biology*, 32(1), 79-82. \doi{10.1007/Bf00160376}
#'
#' **Other References**
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#' *Mathematical Biosciences*, 195(1), 47-55. \doi{10.1016/j.mbs.2005.02.001}
#'
#' @export

SO93F2 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {

  # Determine number of records
  n <- length(records)
  if (init.time == min(records)) {
    n <- n - 1
  }

  # Remove first record, if using as initial time
  if (init.time == min(records)) {
    records <- sort(records)[2:length(records)]
  }

  # Determine length of sighting record
  tn <- max(records) - init.time

  # Determine current time
  bigT <- test.time - init.time

  # Determine s value
  s <- sum(records - init.time)

  # Calculate p-value
  p.value <- Fx(x = tn, s = s, n = n) /
    Fx(x = bigT, s = s, n = n)

  # Calculate numerator for point estimate
  i <- 0:floor(s / tn)
  part1 <- (-1) ^ i
  part2 <- choose(n, i)
  part3 <- (s - (i * tn)) ^ (n - 1)
  numerator <- sum(part1 * part2 * part3)
  rm(i, part1, part2, part3)

  # Calculate denominator for point estimate
  part0 <- n * (n - 1)
  i <- 0:(floor(s / tn) - 1)
  part1 <- (-1) ^ i
  part2 <- choose(n - 1, i)
  part3 <- (s - ((i + 1) * tn)) ^ (n - 2)
  denominator <- part0 * sum(part1 * part2 * part3)
  rm(i, part0, part1, part2, part3)

  # Set up Fopt to help find confidence interval
  Fopt <- function(x) {
    value <- (Fx(x = tn, s = s, n = n) / Fx(x = x, s = s, n = n)) - alpha
    return(value)
  }

  # Numerically find the confidence interval
  conf.int <- uniroot(f = Fopt, interval = c(tn, .Machine$integer.max))$root

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.value = p.value,
    estimate = tn + numerator / denominator,
    conf.int = c(tn, conf.int))

  return(output)

}
