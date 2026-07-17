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
#' @param priors `list` with two elements: `theta` and `lambda`, themselves both
#' `numeric` vectors of length two. The two elements in `theta` are the shape
#' parameters for the Beta hyperprior on \eqn{\theta}. They default to (1, 10).
#' The two elements in `lambda` are the shape and rate parameters for the Gamma
#' prior on \eqn{\lambda}. They default to (1, 1).
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
                   priors = list(theta = c(1, 10), lambda = c(1, 1)),
                   n.chains = 4, n.iter = 11e4, n.burnin = 1e4, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("package 'rjags' is required but could not be found")
  }

  # Check that data and priors are in a valid format
  if (anyNA(records) || any(records < 0) || any(records != floor(records))) {
    stop("records must be non-negative integer counts")
  }
  records <- as.integer(records)

  if (is.null(priors$theta) || length(priors$theta) != 2 ||
    anyNA(priors$theta) || any(priors$theta <= 0)) {
    stop("priors$theta must be a positive vector of length 2")
  }

  if (is.null(priors$lambda) || length(priors$lambda) != 2 ||
    anyNA(priors$lambda) || any(priors$lambda <= 0)) {
    stop("priors$lambda must be a positive vector of length 2")
  }

  # Calculate key values
  bigT <- length(records)
  t_m <- max(which(records > 0))

  # Specify model and parameters
  data_list <- list(
    y = records,
    bigT = bigT,
    t_m = t_m,
    theta_a = priors$theta[1],
    theta_b = priors$theta[2],
    lambda_a = priors$lambda[1],
    lambda_b = priors$lambda[2]
  )

  model_string <- "
    model {
      # 1. Priors
      theta ~ dbeta(theta_a, theta_b)
      tau_L ~ dnegbin(theta, 1)

      lambda ~ dgamma(lambda_a, lambda_b)

      # 2. Likelihood
      for (t in 1:bigT) {
        extant[t] <- step(t_m + tau_L - t)
        mu[t] <- extant[t] * lambda
        y[t] ~ dpois(mu[t])
      }
    }
  "

  inits_list <- function() {
    list(
      theta = runif(1, 0.01, 0.99),
      tau_L = sample(0:(2 * bigT), 1),
      lambda = rgamma(1, priors$lambda[1], priors$lambda[2])
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
      "tau_L", "lambda", "theta"
    ), n.iter = n.iter, thin = n.thin)
  }))

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(samples))
  posterior$time <- init.time + t_m + posterior$tau_L - 1

  # Calculate p(extant)
  p.extant <- mean(posterior$time >= test.time)

  # Calculate point estimate
  estimate <- median(posterior$time)

  # Calculate credible interval bounds
  cred.int.lower <- as.numeric(quantile(posterior$time, 0))
  cred.int.upper <- as.numeric(quantile(posterior$time, 1 - alpha))

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
