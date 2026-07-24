#' @title Backstrom et al.'s (2026) "Certain-only Detectability" model
#'
#' @description
#' Model 3 from Backstrom et al. 2026. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-sided
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `cdis` format (see
#' \code{\link{convert_dodo}} for details).
#' @param detectability a `data.frame` of detectability data, with the same
#' number of rows as `records` and as many columns as there are detection
#' variables.
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the end of the observation period.
#' @param tol tolerance for the lagged extinction hazard tail approximation. The
#' lagged hazard is exactly evaluated until it reaches `tol` of `theta_max`,
#' after which point the remaining tail is approximated using a constant hazard
#' of `theta_max`. Defaults to 0.999; if this yields slow evaluation (in
#' conjunction with a large value of `g`), consider a lower tolerance e.g. 0.99.
#' @param priors `list` with three elements: `theta_max`, `g`, and `sigma`. Both
#' `theta_max` and `sigma` are `numeric` vectors, `g` is a single `numeric.`
#' `theta_max` is of length two, with the two elements being the shape
#' parameters for the Beta hyperprior on \eqn{\theta_{max}}. They default to
#' (1, 10). `g` is the delay parameter for the lagged extinction hazard. It
#' defaults to 0, which recovers the constant-hazard geometric prior for
#' extinction time. `sigma` should either be of length one, or the same length
#' as the number of coefficients to estimate (i.e. `ncol(detectability) + 1`).
#' Defaults to 1.
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

BA26B3 <- function(records, detectability, alpha = 0.05, init.time,
                   test.time = init.time + length(records) - 1, tol = 0.999,
                   priors = list(theta_max = c(1, 10), g = 0, sigma = c(1)),
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

  if (is.null(priors$theta_max) || length(priors$theta_max) != 2 ||
    anyNA(priors$theta_max) || any(priors$theta_max <= 0)) {
    stop("priors$theta_max must be a positive vector of length 2")
  }

  if (is.null(priors$g) || length(priors$g) != 1 ||
    anyNA(priors$g) || priors$g < 0) {
    stop("priors$g must be a single non-negative number")
  }

  detectability <- as.matrix(detectability)
  if (anyNA(detectability) || !is.numeric(detectability)) {
    stop("detectability must be numeric and not contain NA values")
  }

  if (nrow(detectability) != length(detectability)) {
    stop("detectability must have one row per record")
  }

  if (any(is.null(priors$sigma)) || any(priors$sigma <= 0)) {
    stop("priors$sigma must be positive")
  }

  if (length(priors$sigma) != 1 && length(priors$sigma) !=
    (ncol(detectability) + 1)) {
    stop("priors$sigma must be of length 1 or ncol(detectability) + 1")
  }

  # Reformat sigma prior
  if (length(priors$sigma) == 1) {
    sigma <- rep(priors$sigma, ncol(detectability) + 1)
  } else {
    sigma <- priors$sigma
  }

  if (is.null(tol) || !is.numeric(tol) || length(tol) != 1 || is.na(tol) ||
    tol <= 0 || tol >= 1) {
    stop("tol must be a single number in (0, 1)")
  }

  # Calculate key values
  bigT <- length(records)
  t_m <- max(which(records > 0))
  max_lag_obs <- bigT - t_m
  haz_tol <- tol
  if (priors$g == 0) {
    max_lag_prior <- max_lag_obs
  } else {
    lag_req <- ceiling(log(1 - haz_tol) / log(priors$g / (priors$g + 1)) - 1)
    max_lag_prior <- max(max_lag_obs, lag_req)
  }
  n_cat <- max_lag_prior + 2
  y_logfact <- lfactorial(records)
  precision <- 1 / sigma^2

  # Check if the maximum lag before swapping to tail approximation is large
  if (max_lag_prior > 1000) {
    warning(paste0(
      "maximum lag is very large (", max_lag_prior, "); JAGS may be slow"
    ))
  }

  # Specify model and parameters
  data_list <- list(
    y = records,
    x = detectability,
    p = ncol(detectability),
    bigT = bigT,
    t_m = t_m,
    max_lag_prior = max_lag_prior,
    n_cat = n_cat,
    y_logfact = y_logfact,
    zeros = 0L,
    g = priors$g,
    theta_max_a = priors$theta_max[1],
    theta_max_b = priors$theta_max[2],
    precision = precision
  )

  model_string <- "
    model {
      # 1. Priors
      theta_max ~ dbeta(theta_max_a, theta_max_b)

      alpha0 ~ dnorm(0, precision[1])
      for (m in 1:p) {
        alpha[m] ~ dnorm(0, precision[m + 1])
      }

      ## 1.1. Time-varying hazard prior for tau_L up to max_lag_prior
      surv[1] <- 1

      for (k in 1:(max_lag_prior + 1)) {
        lag[k] <- k - 1
        theta_lag[k] <- theta_max * (1 - pow(1 - 1 / (g + 1), lag[k] + 1))
        p_tau[k] <- surv[k] * theta_lag[k]
        surv[k + 1] <- surv[k] * (1 - theta_lag[k])
      }

      ## 1.2. Approximate tail probability after max_lag_prior
      p_tau[n_cat] <- surv[max_lag_prior + 2]
      tau_cat ~ dcat(p_tau[1:n_cat])
      is_tail <- equals(tau_cat, n_cat)

      tau_tail ~ dnegbin(theta_max, 1)

      tau_L <- (1 - is_tail) * (tau_cat - 1) + is_tail *
        (max_lag_prior + 1 + tau_tail)
      tau_E <- t_m + tau_L

      # 2. Likelihood
      for (t in 1:bigT) {
        eta[t] <- alpha0 + inprod(alpha[1:p], x[t, 1:p])
        lambda[t] <- exp(eta[t])

        loglik_obs[t] <- -lambda[t] + y[t] * eta[t] - y_logfact[t]
      }

      cum_loglik[1] <- loglik_obs[1]

      for (t in 2:bigT) {
        cum_loglik[t] <- cum_loglik[t - 1] + loglik_obs[t]
      }

      for (t in 1:bigT) {
        loglik[t] <- cum_loglik[t]
      }

      loglik[bigT + 1] <- cum_loglik[bigT]

      idx <- step(bigT - tau_E) * tau_E + step(tau_E - bigT - 1) * (bigT + 1)
      selected_loglik <- loglik[idx]

      phi <- -selected_loglik
      zeros ~ dpois(phi)
    }
  "

  inits_list <- function() {
    list(
      theta_max = runif(1, 0.01, 0.99),
      tau_cat = sample(1:n_cat, 1),
      tau_tail = sample(0:(2 * bigT), 1),
      alpha0 = rnorm(1, 0, 1),
      alpha = rnorm(ncol(detectability), 0, 1)
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
      "tau_E", "theta_max", "alpha0", "alpha"
    ), n.iter = n.iter, thin = n.thin)
  }))

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(samples))
  posterior$time <- init.time + posterior$tau_E - 1

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
    detectability = detectability,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    tol = tol,
    priors = priors,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(output)
}
