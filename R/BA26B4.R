#' @title Backstrom et al.'s (2026) "Uncertain Effort" model
#'
#' @description
#' Model 4 from Backstrom et al. 2026. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-sided
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `udis` format (see
#' \code{\link{convert_dodo}} for details).
#' @param effort a `data.frame` of effort data, with the same number of rows as
#' `records` and as many columns as there are effort variables.
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the end of the observation period.
#' @param priors `list` with two elements: `sigma_v` and `sigma_i`, the scale
#' parameters for the Normal priors on \emph{\eqn{\beta}} and
#' \emph{\eqn{\gamma}}. Both `sigma_v` and `sigma_i` should either be a single
#' `numeric` object, or a `numeric` vector of the same length as the number of
#' coefficients to estimate (i.e. `ncol(effort) + 1`). Both default to 1.
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
#' This model incorporates both sighting uncertainty and variable survey effort.
#'
#' @references
#' **Key Reference**
#'
#' Backstrom, L. J. et al. (in prep).
#'
#' @seealso [BA26B1()], [BA26B2()], [BA26B3()]
#'
#' @examples
#' \dontrun{
#' # ADD
#' }
#'
#' @export

BA26B4 <- function(records, effort, alpha = 0.05, init.time,
                   test.time = init.time + nrow(records) - 1,
                   priors = list(sigma_v = 1, sigma_i = 1),
                   n.chains = 4, n.iter = 11e4,
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

  effort <- as.matrix(effort)
  if (anyNA(effort) || !is.numeric(effort)) {
    stop("effort must be numeric and not contain NA values")
  }

  if (nrow(effort) != nrow(records)) {
    stop("effort must have one row per record")
  }

  if (any(is.null(priors$sigma_v)) || any(is.null(priors$sigma_i)) ||
    any(priors$sigma_v <= 0) || any(priors$sigma_i <= 0)) {
    stop("both sigma priors must be positive")
  }

  if (length(priors$sigma_v) != 1 && length(priors$sigma_v) !=
    (ncol(effort) + 1) ||
    length(priors$sigma_i) != 1 && length(priors$sigma_i) !=
      (ncol(effort) + 1)) {
    stop("both sigma priors must be of length 1 or ncol(effort) + 1")
  }

  # Reformat sigma priors
  if (length(priors$sigma_v) == 1) {
    sigma_v <- rep(priors$sigma_v, ncol(effort) + 1)
  } else {
    sigma_v <- priors$sigma_v
  }

  if (length(priors$sigma_i) == 1) {
    sigma_i <- rep(priors$sigma_i, ncol(effort) + 1)
  } else {
    sigma_i <- priors$sigma_i
  }

  y_c <- records$certain
  y_u <- records$uncertain

  # Calculate key values
  bigT <- nrow(records)
  p <- ncol(effort)
  logfact_c <- lfactorial(y_c)
  logfact_u <- lfactorial(y_u)
  no_certain_after <- integer(bigT)
  for (t in 1:bigT) {
    if (t == bigT) {
      no_certain_after[t] <- 1L
    } else {
      no_certain_after[t] <- as.integer(all(y_c[(t + 1):bigT] == 0))
    }
  }
  precision_v <- 1 / sigma_v^2
  precision_i <- 1 / sigma_i^2

  # Specify model and parameters
  data_list <- list(
    y_c = y_c,
    y_u = y_u,
    x = effort,
    p = ncol(effort),
    bigT = bigT,
    logfact_c = logfact_c,
    logfact_u = logfact_u,
    no_certain_after = no_certain_after,
    zeros = 0L,
    precision_v = precision_v,
    precision_i = precision_i
  )

  model_string <- "
    model {
      # 1. Priors
      theta ~ dunif(0, 1)
      tau_e ~ dnegbin(theta, 1)
      tau_e1 <- tau_e + 1

      pi_e ~ dunif(0, 1)

      beta0 ~ dnorm(0, precision_v[1])
      gamma0 ~ dnorm(0, precision_i[1])

      for (m in 1:p) {
        beta[m] ~ dnorm(0, precision_v[m + 1])
        gamma[m] ~ dnorm(0, precision_i[m + 1])
      }

      # 2. Likelihood
      for (t in 1:bigT) {
        eta_v[t] <- beta0 + inprod(beta[1:p], x[t, 1:p])
        eta_i[t] <- gamma0 + inprod(gamma[1:p], x[t, 1:p])

        lambda_v[t] <- exp(eta_v[t])
        lambda_i[t] <- exp(eta_i[t])

        mu_c_extant[t] <- lambda_v[t] * pi_e
        mu_u_extant[t] <- lambda_v[t] * (1 - pi_e) + lambda_i[t]

        loglik_c_extant[t] <- -mu_c_extant[t] + y_c[t] *
          (eta_v[t] + log(pi_e)) - logfact_c[t]
        loglik_u_extant[t] <- -mu_u_extant[t] + y_u[t] *
          log(mu_u_extant[t]) - logfact_u[t]

        loglik_extant[t] <- loglik_c_extant[t] + loglik_u_extant[t]

        loglik_u_post[t] <- -lambda_i[t] +  y_u[t] * eta_i[t] - logfact_u[t]
      }

      cum_extant[1] <- loglik_extant[1]
      cum_u_post[1] <- loglik_u_post[1]

      for (t in 2:bigT) {
        cum_extant[t] <- cum_extant[t - 1] + loglik_extant[t]
        cum_u_post[t] <- cum_u_post[t - 1] + loglik_u_post[t]
      }

      for (t in 1:bigT) {
        post_u_after[t] <- cum_u_post[bigT] - cum_u_post[t]

        loglik_raw[t] <- cum_extant[t] + post_u_after[t]

        loglik[t] <- no_certain_after[t] * loglik_raw[t] +
          (1 - no_certain_after[t]) * (-1.0E12)
      }
      loglik[bigT + 1] <- cum_extant[bigT]

      idx <- step(bigT - tau_e1) * tau_e1 + step(tau_e1 - bigT - 1) * (bigT + 1)

      phi <- -loglik[idx]
      zeros ~ dpois(phi)
    }
  "

  inits_list <- function() {
    list(
      theta = runif(1, 0.01, 0.99),
      tau_e = sample(0:(2 * bigT), 1),
      pi_e = runif(1, 0.01, 0.99),
      beta0 = rnorm(1, 0, 1),
      gamma0 = rnorm(1, 0, 1),
      beta = rnorm(ncol(effort), 0, 1),
      gamma = rnorm(ncol(effort), 0, 1)
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
      "tau_e1", "pi_e", "theta", "beta0", "beta", "gamma0", "gamma"
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
