#' @title Weiss & Marshall's (1999) "Discrete-time" model
#'
#' @description
#' The model from Weiss & Marshall 1999. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-tailed
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `cbin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param surveys a `numeric` vector of the survey times for each observation
#' in `records`.
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at.
#' @param priors `list` with three elements: `lambda`, `c` and `d`. `lambda` is
#' the mean lifetime (half-life) for the exponential prior on `S`, the time of
#' extinction. `c` and `d` are the two shape parameters for the beta prior on
#' `pi`, the pre-extinction detection probability.
#'
#' @returns a `list` object with the original parameters and the p(extant),
#' point estimate, and credible interval included as elements. The credible
#' interval is a two-element numeric vector called `cred.int`. The point
#' estimate is the median (not the mean) of the posterior distribution of
#' extinction time.
#'
#' @note
#' All sighting records are assumed to be certain.
#'
#' @references
#' **Key Reference**
#'
#' Weiss, R. E., & Marshall, C. R. (1999). The Uncertainty in the True End
#' Point of a Fossil's Stratigraphic Range When Stratigraphic Sections Are
#' Sampled Discretely. *Mathematical Geology*, 31(4), 435-453.
#' \doi{10.1023/A:1007542725180}
#'
#' @examples
#' # Run the Verneuilinoides sp. A analysis from Weiss & Marshall 1999
#' WM99B1(verneuilinoides, weissmarshall_surveys,
#'   priors = list(lambda = 800, c = 11, d = 70)
#' )
#' # Run the Eggerellina brevis analysis from Weiss & Marshall 1999
#' WM99B1(eggerellina_brevis, weissmarshall_surveys,
#'   priors = list(lambda = 800, c = 85, d = 17)
#' )
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' WM99B1(curlew, 1817:2022, priors = list(lambda = 1e6, c = 1, d = 1))
#' }
#'
#' @export

WM99B1 <- function(records, surveys, alpha = 0.05, init.time = min(surveys),
                   test.time = max(surveys), priors) {
  # Sort records and surveys
  rs <- data.frame(
    t = surveys,
    y = records
  )

  rs <- sort_by(rs, ~t)
  t <- rs$t
  y <- rs$y
  rm(rs)

  t0 <- min(t)
  t <- t - t0

  # Calculate key values
  m <- max(which(y == 1))
  tm <- t[m]
  am <- sum(y)
  J <- length(t)

  # Helper functions
  check <- function(S, tm) {
    return(S >= tm)
  }

  aS <- function(S, y, t) {
    return(sum(y[which(t <= S)]))
  }

  bS <- function(S, y, t) {
    return(sum(1 - y[which(t <= S)]))
  }

  gj <- function(am, bj, c, d, lambda, tj, tj1) {
    gj <- exp(lgamma(am + c) + lgamma(bj + d) - lgamma(am + bj + c + d) +
      log(exp(-lambda^(-1) * tj) - exp(-lambda^(-1) * tj1)))
    return(gj)
  }

  ding <- function(a, t0, t1) {
    f <- function(x, a) {
      return(x^(a - 1) * exp(-x))
    }
    integrand <- integrate(f, t0, t1, a = a)
    return(gamma(a)^(-1) * integrand$value)
  }

  # Calculate constant of proportionality
  K0 <- c()
  for (j in m:J) {
    K0[j] <- gj(
      am, bS(t[j], y, t), priors$c, priors$d, priors$lambda, t[j],
      ifelse(is.na(t[j + 1]), Inf, t[j + 1])
    )
  }
  K0 <- sum(K0, na.rm = TRUE)

  # Posterior quantiles
  q <- c()
  for (l in m:J) {
    sums <- c()
    for (j in m:l) {
      sums[j] <- K0^(-1) * gj(
        am, bS(t[j], y, t), priors$c, priors$d,
        priors$lambda, t[j], ifelse(is.na(t[j + 1]),
          Inf, t[j + 1]
        )
      )
    }
    q[l + 1] <- sum(sums, na.rm = TRUE)
  }

  pq <- function(t, q, epsilon, lambda) {
    l <- min(which(q > epsilon)) - 1
    tl <- t[l]
    tl1 <- ifelse(is.na(t[l + 1]), Inf, t[l + 1])
    ql <- ifelse(is.na(q[l]), 0, q[l])
    ql1 <- q[l + 1]

    return(-lambda * log(exp(-tl / lambda) - (
      exp(-tl / lambda) - exp(-tl1 / lambda)) * ((epsilon - ql) / (ql1 - ql))))
  }

  # Get quantiles
  q050 <- pq(t, q, 0.5, priors$lambda)
  q950 <- pq(t, q, 0.95, priors$lambda)

  # Get p(extinct)
  pqe <- function(t, q, epsilon, lambda) {
    return(pq(t, q, epsilon, lambda) - test.time + t0)
  }

  p.extinct <- uniroot(
    f = pqe, interval = c(0, 1 - 1e-9), t = t, q = q,
    lambda = priors$lambda
  )$root

  # Output
  output <- list(
    records = records,
    surveys = surveys,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    priors = priors,
    p.extant = 1 - p.extinct,
    estimate = q050 + t0,
    cred.int = c(tm + t0, q950 + t0)
  )

  return(output)
}
