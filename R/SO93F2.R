#' @title Solow's (1993) "Declining" model
#'
#' @description
#' Equations 3 and 4 from Solow 1993 (Equation 7 from Solow 2005), and
#' Equation 8 from Solow 2005. Estimates a p-value for testing competing
#' hypotheses of extinction/non-extinction, and a one-sided \eqn{1 - \alpha}
#' confidence interval and point estimate on the time of extinction.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case  this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param precBits number of bits of precision to use in `Rmpfr` arithmetic.
#' Defaults to 1024, which should be sufficient for most datasets.
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
#' @examples
#' # Run the Black-footed Ferret analysis from Solow 1993
#' SO93F2(ferret$ccon, init.time = 0, test.time = 229)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' SO93F2(curlew$ccon, init.time = 1817, test.time = 2022)
#' }
#'
#' @export

SO93F2 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   precBits = 1024) {
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

  # Determine s value
  s <- sum(records - init.time)

  # Convert everything to mpfr
  s <- Rmpfr::mpfr(s, precBits = precBits)
  n <- Rmpfr::mpfr(n, precBits = precBits)
  tn <- Rmpfr::mpfr(tn, precBits = precBits)
  bigT <- Rmpfr::mpfr(bigT, precBits = precBits)

  # Calculate p-value
  p.value <- as.numeric(Fx(x = tn, s = s, n = n) /
    Fx(x = bigT, s = s, n = n))

  # Calculate numerator for point estimate
  i <- 0:as.integer(floor(s / tn))
  i <- Rmpfr::mpfr(i, precBits = precBits)
  part1 <- (Rmpfr::mpfr(-1, precBits = precBits))^i
  part2 <- Rmpfr::chooseMpfr(n, i)
  part3 <- (s - i * tn)^(n - 1)
  numerator <- sum(part1 * part2 * part3)
  rm(i, part1, part2, part3)

  # Calculate denominator for point estimate
  part0 <- n * (n - 1)
  i <- 0:as.integer((floor(s / tn) - 1))
  i <- Rmpfr::mpfr(i, precBits = precBits)
  part1 <- (Rmpfr::mpfr(-1, precBits = precBits))^i
  part2 <- Rmpfr::chooseMpfr(n - 1, i)
  part3 <- (s - (i + 1) * tn)^(n - 2)
  denominator <- part0 * sum(part1 * part2 * part3)
  rm(i, part0, part1, part2, part3)

  # Set up Fopt to help find confidence interval
  Fopt <- function(x) {
    value <- as.numeric(Fx(x = tn, s = s, n = n) / Fx(x = x, s = s, n = n)) -
      alpha
    return(value)
  }

  # Numerically find the confidence interval
  conf.int <- tryCatch(
    uniroot(
      f = Fopt,
      interval = c(as.integer(tn), as.integer(s))
    )$root,
    error = function(e) NA
  )

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    precBits = precBits,
    p.value = p.value,
    estimate = as.numeric(init.time + tn + numerator / denominator),
    conf.int = c(init.time + as.numeric(tn), init.time + conf.int)
  )

  return(output)
}

#' @title F(x) from Solow 2005
#'
#' @description
#' Helper function. Modified Equation 7 from Solow 2005.
#'
#' @param x generally \eqn{t_n} or \eqn{T}.
#' @param s sum of relative sighting times.
#' @param n number of sighting records.
#'
#' @returns a number.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#' *Mathematical Biosciences*, 195(1), 47-55. \doi{10.1016/j.mbs.2005.02.001}
#'
#' @noRd

Fx <- function(x, s, n) {
  if (as.integer(floor(s / x)) == 0) {
    warning("floor(s / x) is zero - NA produced")
    return(NA)
  }

  y <- x / s
  is <- 1:as.integer(floor(1 / y))

  vals <- list()
  for (i in is) {
    vals[[i]] <- (-1)^(i - 1) * Rmpfr::chooseMpfr(n, i) * (1 - i * y)^(n - 1)
  }

  return(1 - sum(do.call(c, vals)))
}
