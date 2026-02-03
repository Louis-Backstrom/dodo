#' @title Weiss & Marshall's (1999) "Discrete-time" model
#'
#' @description
#' The model from Weiss & Marshall 1999. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-sided
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `cbin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param surveys a `numeric` vector of the survey times for each observation
#' in `records`.
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the time of the final survey.
#' @param priors `list` with three elements: `lambda`, `c` and `d`. `lambda` is
#' the mean lifetime (half-life) for the exponential prior on `S`, the time of
#' extinction. `c` and `d` are the two shape parameters for the beta prior on
#' `pi`, the pre-extinction detection probability.
#' @param increment step size used for integration. Defaults to 0.001.
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
#' # 737.1544 - 7.4 = 729.7544 ≈ 730. from paper
#' # Run the Eggerellina brevis analysis from Weiss & Marshall 1999
#' WM99B1(eggerellina_brevis, weissmarshall_surveys,
#'   priors = list(lambda = 800, c = 85, d = 17)
#' )
#' # 12.18193 - 11.7 = 0.48193 ≈ .482 from paper
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' WM99B1(curlew$cbin, 1817:2022, priors = list(lambda = 1e3, c = 1, d = 1),
#'  increment = 0.01)
#' }
#'
#' @export

WM99B1 <- function(records, surveys, alpha = 0.05, test.time = max(surveys),
                   priors, increment = 0.001) {
  # Sort records and surveys
  rs <- data.frame(
    t = surveys,
    y = records
  )

  rs <- sort_by(rs, ~t)
  t <- rs$t
  y <- rs$y

  t0 <- min(t)
  t <- t - t0

  # Calculate key values
  m <- max(which(y == 1))
  tm <- t[m]
  am <- sum(y)
  J <- length(t)

  # Helper functions
  bS <- function(S, y, t) {
    return(sum(1 - y[t <= S]))
  }

  gj <- function(am, bj, c, d, lambda, tj, tj1) {
    gj <- exp(lgamma(am + c) + lgamma(bj + d) - lgamma(am + bj + c + d) +
                log(exp(-lambda^(-1) * tj) - exp(-lambda^(-1) * tj1)))
    return(gj)
  }

  # Precompute cumulative failures
  cum_surveys <- seq_len(J)
  cum_detections <- cumsum(y)
  cum_failures <- cum_surveys - cum_detections

  # Calculate constant of proportionality
  K0 <- numeric(J)
  for (j in m:J) {
    K0[j] <- gj(
      am, cum_failures[j], priors$c, priors$d, priors$lambda, t[j],
      ifelse(is.na(t[j + 1]), Inf, t[j + 1])
    )
  }
  K0_sum <- sum(K0, na.rm = TRUE)

  # Posterior quantiles
  q <- numeric(J + 1)
  for (l in m:J) {
    q[l + 1] <- sum(K0[m:l]) / K0_sum
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
  qupper <- pq(t, q, 1 - alpha, priors$lambda)

  # Get p(extinct)
  pqe <- function(t, q, epsilon, lambda) {
    return(pq(t, q, epsilon, lambda) - test.time + t0)
  }

  # Compute pq_min and pq_max
  pq_min <- pq(t, q, 0, priors$lambda)
  pq_max <- pq(t, q, 1 - 1e-9, priors$lambda)

  # Check edge cases and calculate p(extinct)
  if (test.time > pq_max) {
    # Too close to 1
    p.extinct <- 1
  } else if (test.time < pq_min) {
    # Too close to 0
    p.extinct <- 0
  } else {
    # In meaningful [0, 1] range
    p.extinct <- uniroot(
      f = pqe, interval = c(0, 1 - 1e-9), t = t, q = q,
      lambda = priors$lambda
    )$root
  }

  # Get posterior mean
  S <- seq(tm, t[J] + 10 * priors$lambda, by = increment) # estimate
  bS_S <- cum_failures[pmin(findInterval(S, t), J)]

  S_posterior <- exp(
    -log(K0_sum) +
      lgamma(am + priors$c) +
      lgamma(bS_S + priors$d) -
      lgamma(am + bS_S + priors$c + priors$d) +
      log(priors$lambda^(-1)) +
      (-priors$lambda^(-1) * S)
  )

  mean <- sum(S * S_posterior) / sum(S_posterior)

  # Output
  output <- list(
    records = records,
    surveys = surveys,
    alpha = alpha,
    test.time = test.time,
    priors = priors,
    increment = increment,
    p.extant = 1 - p.extinct,
    estimate = mean + t0,
    cred.int = c(tm + t0, qupper + t0)
  )

  return(output)
}


