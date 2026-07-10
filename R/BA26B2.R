#' @title Backstrom et al.'s (2026) "Uncertain No-effort" model
#'
#' @description
#' Model 2 from Backstrom et al. 2026. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-sided
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `udis` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the end of the observation period.
#' @param priors `list` with two elements: `a` and `b`, the shape and rate
#' parameters for the Gamma prior on \eqn{\lambda}. Both default to 1.
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
#' Sampling effort is assumed to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Backstrom, L. J. et al. (in prep).
#'
#' @seealso [BA26B1()], [BA26B3()], [BA26B4()]
#'
#' @examples
#' \dontrun{
#' # ADD
#' }
#'
#' @export

BA26B2 <- function(records, alpha = 0.05, init.time,
                   test.time = init.time + nrow(records) - 1,
                   priors = list(a = 1, b = 1), n.chains = 4, n.iter = 11e4,
                   n.burnin = 1e4, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("package 'rjags' is required but could not be found")
  }

  # Check that data and priors are in a valid format
  if (anyNA(records$certain) || anyNA(records$uncertain) ||
    any(records$certain < 0) || any(records$uncertain < 0) ||
    any(records$certain != floor(records$certain)) ||
    any(records$uncertain != floor(records$uncertain))) {
    stop("records must be non-negative integer counts")
  }
  records$certain <- as.integer(records$certain)
  records$uncertain <- as.integer(records$uncertain)

  if (is.null(priors$a) || is.null(priors$b) ||
    priors$a <= 0 || priors$b <= 0) {
    stop("priors$a and priors$b must be positive")
  }

  y_c <- records$certain
  y_u <- records$uncertain

  # Calculate key values
  bigT <- nrow(records)
  y_c_sum <- cumsum(y_c)
  y_c_logfact_sum <- cumsum(lfactorial(y_c))
  y_c_total <- sum(y_c)
  y_c_logfact_total <- sum(lfactorial(y_c))
  y_u_sum <- cumsum(y_u)
  y_u_logfact_sum <- cumsum(lfactorial(y_u))
  y_u_total <- sum(y_u)
  y_u_logfact_total <- sum(lfactorial(y_u))
  no_certain_after <- integer(bigT)
  for (t in 1:bigT) {
    if (t == bigT) {
      no_certain_after[t] <- 1L
    } else {
      no_certain_after[t] <- as.integer(all(y_c[(t + 1):bigT] == 0))
    }
  }

  # Specify model and parameters
  data_list <- list(
    bigT = bigT,
    y_c_sum = y_c_sum,
    y_c_logfact_sum = y_c_logfact_sum,
    y_c_total = y_c_total,
    y_c_logfact_total = y_c_logfact_total,
    y_u_sum = y_u_sum,
    y_u_logfact_sum = y_u_logfact_sum,
    y_u_total = y_u_total,
    y_u_logfact_total = y_u_logfact_total,
    no_certain_after = no_certain_after,
    zeros = 0L,
    a = priors$a,
    b = priors$b
  )

  model_string <- "
    model {
      # 1. Priors
      theta ~ dbeta(0.5, 0.5) # Jeffrey's prior
      tau_e ~ dnegbin(theta, 1)
      tau_e1 <- tau_e + 1

      lambda_v ~ dgamma(a, b)
      lambda_i ~ dgamma(a, b)

      pi_e ~ dunif(0, 1)

      # 2. Likelihood
      mu_c_extant <- lambda_v * pi_e
      mu_u_extant <- lambda_v * (1 - pi_e) + lambda_i

      for (t in 1:bigT) {
        n_after[t] <- bigT - t
        y_u_after[t] <- y_u_total - y_u_sum[t]
        y_u_logfact_after[t] <- y_u_logfact_total - y_u_logfact_sum[t]

        loglik_c_extant[t] <- -t * mu_c_extant + y_c_sum[t] *
          log(mu_c_extant) - y_c_logfact_sum[t]
        loglik_u_extant[t] <- -t * mu_u_extant + y_u_sum[t] *
          log(mu_u_extant) - y_u_logfact_sum[t]
        loglik_u_post[t] <- -n_after[t] * lambda_i + y_u_after[t] *
          log(lambda_i) - y_u_logfact_after[t]

        loglik_raw[t] <- loglik_c_extant[t] + loglik_u_extant[t] +
          loglik_u_post[t]
        loglik[t] <- no_certain_after[t] * loglik_raw[t] +
          (1 - no_certain_after[t]) * (-1.0E12)
      }

      mu_c_extant_after <- lambda_v * pi_e
      mu_u_extant_after <- lambda_v * (1 - pi_e) + lambda_i

      loglik[bigT + 1] <- -bigT * mu_c_extant_after + y_c_total *
        log(mu_c_extant_after) - y_c_logfact_total + -bigT * mu_u_extant_after +
        y_u_total * log(mu_u_extant_after) - y_u_logfact_total

      x <- step(bigT - tau_e1) * tau_e1 + step(tau_e1 - bigT - 1) * (bigT + 1)

      phi <- -loglik[x]
      zeros ~ dpois(phi)
    }
  "

  inits_list <- function() {
    list(
      theta = runif(1, 0.01, 0.99),
      tau_e = sample(0:(2 * bigT), 1),
      lambda_v = rgamma(1, shape = priors$a, rate = priors$b),
      lambda_i = rgamma(1, shape = priors$a, rate = priors$b),
      pi_e = runif(1, 0.01, 0.99)
    )
  }

  model_file <- tempfile(fileext = ".txt")
  writeLines(model_string, con = model_file)
  on.exit(unlink(model_file), add = TRUE)

  # Run MCMC sampling
  invisible(capture.output({
    jags_model <- rjags::jags.model(
      file = model_file, data = data_list, inits = inits_list,
      n.chains = n.chains, n.adapt = n.burnin
    )
    update(jags_model, n.iter = n.burnin)
    samples <- rjags::coda.samples(jags_model, variable.names = c(
      "tau_e1", "lambda_v", "lambda_i", "pi_e", "theta"
    ), n.iter = n.iter, thin = n.thin)
  }))

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(samples))
  posterior$year <- posterior$tau_e1 + init.time - 1

  # Calculate p(extant)
  p.extant <- mean(posterior$year >= test.time)

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
    priors = priors,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(output)
}
