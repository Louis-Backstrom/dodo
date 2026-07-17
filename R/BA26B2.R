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
#' @param priors `list` with three elements: `theta`, `lambda_v` and `lambda_i`,
#' themselves all `numeric` vectors of length two. The two elements in `theta`
#' are the shape parameters for the Beta hyperprior on \eqn{\theta}. They
#' default to (1, 10). The two elements in `lambda_v` and `lambda_i` are the
#' shape and rate parameters for the Gamma priors on \eqn{\lambda_v} and
#' \eqn{\lambda_i}. They both default to (1, 1).
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
                   priors = list(
                     theta = c(1, 10), lambda_v = c(1, 1), lambda_i = c(1, 1)
                   ),
                   n.chains = 4, n.iter = 11e4, n.burnin = 1e4, n.thin = 10) {
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

  if (is.null(priors$theta) || length(priors$theta) != 2 ||
    anyNA(priors$theta) || any(priors$theta <= 0)) {
    stop("priors$theta must be a positive vector of length 2")
  }

  if (is.null(priors$lambda_v) || length(priors$lambda_v) != 2 ||
    anyNA(priors$lambda_v) || any(priors$lambda_v <= 0)) {
    stop("priors$lambda_v must be a positive vector of length 2")
  }

  if (is.null(priors$lambda_i) || length(priors$lambda_i) != 2 ||
    anyNA(priors$lambda_i) || any(priors$lambda_i <= 0)) {
    stop("priors$lambda_i must be a positive vector of length 2")
  }

  y_c <- records$certain
  y_u <- records$uncertain

  # Calculate key values
  bigT <- nrow(records)
  t_m <- max(which(y_c > 0))


  # Specify model and parameters
  data_list <- list(
    y_c = y_c,
    y_u = y_u,
    bigT = bigT,
    t_m = t_m,
    theta_a = priors$theta[1],
    theta_b = priors$theta[2],
    lambda_v_a = priors$lambda_v[1],
    lambda_v_b = priors$lambda_v[2],
    lambda_i_a = priors$lambda_i[1],
    lambda_i_b = priors$lambda_i[2]
  )

  model_string <- "
    model {
      # 1. Priors
      theta ~ dbeta(theta_a, theta_b)
      tau_L ~ dnegbin(theta, 1)

      lambda_v ~ dgamma(lambda_v_a, lambda_v_b)
      lambda_i ~ dgamma(lambda_i_a, lambda_i_b)

      pi_e ~ dunif(0, 1)

      # 2. Likelihood
      for (t in 1:bigT) {
        extant[t] <- step(t_m + tau_L - t)

        # Certain
        mu_c[t] <- extant[t] * lambda_v * pi_e
        y_c[t] ~ dpois(mu_c[t])

        # Uncertain
        mu_u[t] <- extant[t] * lambda_v * (1 - pi_e) + lambda_i
        y_u[t] ~ dpois(mu_u[t])
      }
    }
  "

  inits_list <- function() {
    list(
      theta = runif(1, 0.01, 0.99),
      tau_L = sample(0:(2 * bigT), 1),
      lambda_v = rgamma(1, priors$lambda_v[1], priors$lambda_i[2]),
      lambda_i = rgamma(1, priors$lambda_v[1], priors$lambda_i[2]),
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
      "tau_L", "lambda_v", "lambda_i", "pi_e", "theta"
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
