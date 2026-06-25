#' @title Backstrom et al.'s (2026) "Certain-only No-effort" model
#'
#' @description
#' Model 1 from Backstrom et al. 2026. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-sided
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `cdis` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the end of the observation period.
#' @param n.chains number of MCMC chains to run. Defaults to 4.
#' @param n.iter number of iterations in each chain. Defaults to 110,000.
#' @param n.burnin number of iterations to discard as burn-in. Defaults to
#' 10,000.
#' @param n.thin thinning rate. Defaults to 10.
#'
#' @returns a `list` object with the original parameters and the p(extant),
#' point estimate, and credible interval included as elements. The credible
#' interval is a two-element numeric vector called `cred.int`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Backstrom, L. J. et al. (in prep).
#'
#' @seealso [BA26B2()], [BA26B3()], [BA26B4()]
#'
#' @examples
#' \dontrun{
#' # ADD
#' }
#'
#' @export

BA26B1 <- function(records, alpha = 0.05, init.time,
                   test.time = init.time + length(records) - 1,
                   n.chains = 4, n.iter = 11e4, n.burnin = 1e4, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("package 'rjags' is required but could not be found")
  }

  t_m <- max(which(records > 0))
  bigT <- length(records)
  y_sum <- cumsum(records)
  lfact_sum <- cumsum(lfactorial(records))
  zero_ok <- integer(bigT)
  for (t in 1:bigT) {
    zero_ok[t] <- as.integer(if (t == bigT) TRUE else all(records[(t + 1):bigT] == 0))
  }

  a <- 1
  b <- 1

  data_list <- list(
    bigT = bigT,
    y_sum = y_sum,
    lfact_sum = lfact_sum,
    zero_ok = zero_ok,
    zeros = 0L,
    a = a,
    b = b
  )

  model_string <- "
    model {
      # 1. Priors
      theta ~ dunif(0, 1)
      tau_e ~ dnegbin(theta, 1)
      tau_e1 <- tau_e + 1

      lambda ~ dgamma(a, b)

      # 2. Likelihood
      for (t in 1:bigT) {
        loglik_raw[t] <- -t * lambda + y_sum[t] * log(lambda) - lfact_sum[t]
        loglik[t] <- zero_ok[t] * loglik_raw[t] + (1 - zero_ok[t]) * (-1.0E12)
      }

      loglik[bigT + 1] <- -bigT * lambda + y_sum[bigT] * log(lambda) -
        lfact_sum[bigT]
      x <- step(bigT - tau_e1) * tau_e1 + step(tau_e1 - bigT - 1) * (bigT + 1)

      phi <- -loglik[x]
      zeros ~ dpois(phi)
    }
  "

  inits_list <- function() {
    list(
      theta = runif(1, 0.01, 0.99),
      tau_e = sample(0:(2 * bigT), 1),
      lambda = rgamma(1, shape = a, rate = b)
    )
  }

  model_file <- tempfile(fileext = ".txt")
  writeLines(model_string, con = model_file)

  invisible(capture.output({
    jags_model <- rjags::jags.model(
      file = model_file, data = data_list, inits = inits_list,
      n.chains = n.chains, n.adapt = n.burnin
    )
    update(jags_model, n.iter = n.burnin)
    samples <- rjags::coda.samples(jags_model, variable.names = c(
      "tau_e1", "lambda", "theta"
    ), n.iter = n.iter, thin = n.thin)
  }))

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(samples))
  posterior$year <- posterior$tau_e1 + init.time - 1

  # Calculate p(extant)
  p.extant <- mean(posterior$year > test.time)

  # Calculate point estimate
  estimate <- median(posterior$year)

  # Calculate credible interval bounds
  cred.int.lower <- as.numeric(quantile(posterior$year, 0))
  cred.int.upper <- as.numeric(quantile(posterior$year, 1 - alpha))

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(output)
}
