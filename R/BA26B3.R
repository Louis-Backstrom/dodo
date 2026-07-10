#' @title Backstrom et al.'s (2026) "Certain-only Effort" model
#'
#' @description
#' Model 3 from Backstrom et al. 2026. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-sided
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `cdis` format (see
#' \code{\link{convert_dodo}} for details).
#' @param effort a `data.frame` of effort data, with the same number of rows as
#' `records` and as many columns as there are effort variables.
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the end of the observation period.
#' @param priors `list` with one element: `sigma`, the scale parameter(s) for
#' the Normal priors on \emph{\eqn{\alpha}}. `sigma` should either be a single
#' `numeric` object, or a `numeric` vector of the same length as the number of
#' coefficients to estimate (i.e. `ncol(effort) + 1`). Defaults to 1.
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
#' All sighting records are assumed to be certain.
#'
#' @references
#' **Key Reference**
#'
#' Backstrom, L. J. et al. (in prep).
#'
#' @seealso [BA26B1()], [BA26B2()], [BA26B4()]
#'
#' @examples
#' \dontrun{
#' # ADD
#' }
#'
#' @export

BA26B3 <- function(records, effort, alpha = 0.05, init.time,
                   test.time = init.time + length(records) - 1,
                   priors = list(sigma = 1), n.chains = 4, n.iter = 11e4,
                   n.burnin = 1e4, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("package 'rjags' is required but could not be found")
  }

  # Check that data and priors are in a valid format
  if (anyNA(records) || any(records < 0) || any(records != floor(records))) {
    stop("records must be non-negative integer counts")
  }
  records <- as.integer(records)

  effort <- as.matrix(effort)
  if (anyNA(effort) || !is.numeric(effort)) {
    stop("effort must be numeric and not contain NA values")
  }

  if (nrow(effort) != length(records)) {
    stop("effort must have one row per record")
  }

  if (any(is.null(priors$sigma)) || any(priors$sigma <= 0)) {
    stop("priors$sigma must be positive")
  }

  if (length(priors$sigma) != 1 && length(priors$sigma) != (ncol(effort) + 1)) {
    stop("priors$sigma must be of length 1 or ncol(effort) + 1")
  }

  # Reformat sigma prior
  if (length(priors$sigma) == 1) {
    sigma <- rep(priors$sigma, ncol(effort) + 1)
  } else {
    sigma <- priors$sigma
  }

  # Calculate key values
  bigT <- length(records)
  logfact_y <- lfactorial(records)
  no_records_after <- integer(bigT)
  for (t in 1:bigT) {
    if (t == bigT) {
      no_records_after[t] <- 1L
    } else {
      no_records_after[t] <- as.integer(all(records[(t + 1):bigT] == 0))
    }
  }
  precision <- 1 / sigma^2

  # Specify model and parameters
  data_list <- list(
    y = records,
    x = effort,
    p = ncol(effort),
    bigT = bigT,
    logfact_y = logfact_y,
    no_records_after = no_records_after,
    zeros = 0L,
    precision = precision
  )

  model_string <- "
    model {
      # 1. Priors
      theta ~ dbeta(0.5, 0.5) # Jeffrey's prior
      tau_e ~ dnegbin(theta, 1)
      tau_e1 <- tau_e + 1

      alpha0 ~ dnorm(0, precision[1])
      for (m in 1:p) {
        alpha[m] ~ dnorm(0, precision[m + 1])
      }

      # 2. Likelihood
      for (t in 1:bigT) {
        eta[t] <- alpha0 + inprod(alpha[1:p], x[t, 1:p])
        lambda[t] <- exp(eta[t])

        loglik_obs[t] <- -lambda[t] + y[t] * eta[t] - logfact_y[t]
      }

      cum_loglik[1] <- loglik_obs[1]
      for (t in 2:bigT) {
        cum_loglik[t] <- cum_loglik[t - 1] + loglik_obs[t]
      }

      for (t in 1:bigT) {
        loglik[t] <- no_records_after[t] * cum_loglik[t] +
          (1 - no_records_after[t]) * (-1.0E12)
      }
      loglik[bigT + 1] <- cum_loglik[bigT]

      idx <- step(bigT - tau_e1) * tau_e1 + step(tau_e1 - bigT - 1) * (bigT + 1)

      phi <- -loglik[idx]
      zeros ~ dpois(phi)
    }
  "

  inits_list <- function() {
    list(
      theta = runif(1, 0.01, 0.99),
      tau_e = sample(0:(2 * bigT), 1),
      alpha0 = rnorm(1, 0, 1),
      alpha = rnorm(ncol(effort), 0, 1)
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
      "tau_e1", "theta", "alpha0", "alpha"
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
    effort = effort,
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
