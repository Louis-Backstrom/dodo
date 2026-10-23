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
#' @param tol tolerance for the lagged extinction hazard tail approximation. The
#' lagged hazard is exactly evaluated until it reaches `tol` of `theta_max`,
#' after which point the remaining tail is approximated using a constant hazard
#' of `theta_max`. Defaults to 0.999; if this yields slow evaluation (in
#' conjunction with a large value of `g`), consider a lower tolerance e.g. 0.99.
#' @param priors `list` with four elements: `theta_max`, `g`, `lambda_v`, and
#' `lambda_i`. `theta_max`, `lambda_v`, and `lambda_i` are all `numeric` vectors
#' of length two, `g` is a single `numeric`. The two elements in `theta_max`
#' are the shape parameters for the Beta hyperprior on \eqn{\theta_{max}}. They
#' default to (1, 10). `g` is the delay parameter for the lagged extinction
#' hazard. It defaults to 0, which recovers the constant-hazard geometric prior
#' for extinction time. The two elements in `lambda_v` and `lambda_i` are the
#' shape and rate parameters for the Gamma priors on \eqn{\lambda_v} and
#' \eqn{\lambda_i}. They both default to (1, 1).
#' @param n.chains number of MCMC chains to run. Defaults to 4.
#' @param n.iter number of iterations in each chain. Defaults to 110,000.
#' @param n.burnin number of iterations to discard as burn-in. Defaults to
#' 10,000.
#' @param n.thin thinning rate. Defaults to 10.
#'
#' @returns a `dodo_fit` object with the original data, parameters, posterior
#' samples, and the p(extant), point estimate, and credible interval included as
#' elements.
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
                   test.time = init.time + nrow(records) - 1, tol = 0.999,
                   priors = list(
                     theta_max = c(1, 10), g = 0,
                     lambda_v = c(1, 1), lambda_i = c(1, 1)
                   ),
                   n.chains = 4, n.iter = 11e4, n.burnin = 1e4, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("package 'rjags' is required but could not be found")
  }

  # Check that scalar arguments are in a valid format
  if (!is.numeric(alpha) || length(alpha) != 1L ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1)")
  }

  if (is.null(tol) || !is.numeric(tol) || length(tol) != 1L ||
    !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("tol must be a single number in (0, 1)")
  }

  # Check that records are in a valid format
  if (!is.data.frame(records) ||
    !all(c("certain", "uncertain") %in% names(records))) {
    stop("records must contain columns named 'certain' and 'uncertain'")
  }

  if (!is.numeric(records$certain) || !is.numeric(records$uncertain) ||
    anyNA(records$certain) || anyNA(records$uncertain) ||
    any(!is.finite(records$certain)) ||
    any(!is.finite(records$uncertain)) ||
    any(records$certain < 0) || any(records$uncertain < 0) ||
    any(records$certain != floor(records$certain)) ||
    any(records$uncertain != floor(records$uncertain))) {
    stop("records must contain only finite, non-negative integer counts")
  }

  if (nrow(records) < 2L) {
    stop("records must contain at least two sampling intervals")
  }

  if (!any(records$certain > 0)) {
    stop("records must contain at least one positive certain count")
  }

  records$certain <- as.integer(records$certain)
  records$uncertain <- as.integer(records$uncertain)

  # Check that extinction time priors are in a valid format
  if (is.null(priors$theta_max) || !is.numeric(priors$theta_max) ||
    length(priors$theta_max) != 2L || anyNA(priors$theta_max) ||
    any(!is.finite(priors$theta_max)) || any(priors$theta_max <= 0)) {
    stop("priors$theta_max must be a finite positive vector of length 2")
  }

  if (is.null(priors$g) || !is.numeric(priors$g) ||
    length(priors$g) != 1L || !is.finite(priors$g) ||
    priors$g < 0) {
    stop("priors$g must be a single finite non-negative number")
  }

  # Check that sighting-rate priors are in a valid format
  if (is.null(priors$lambda_v) || !is.numeric(priors$lambda_v) ||
    length(priors$lambda_v) != 2L || anyNA(priors$lambda_v) ||
    any(!is.finite(priors$lambda_v)) ||
    any(priors$lambda_v <= 0)) {
    stop("priors$lambda_v must be a finite positive vector of length 2")
  }

  if (is.null(priors$lambda_i) || !is.numeric(priors$lambda_i) ||
    length(priors$lambda_i) != 2L || anyNA(priors$lambda_i) ||
    any(!is.finite(priors$lambda_i)) ||
    any(priors$lambda_i <= 0)) {
    stop("priors$lambda_i must be a finite positive vector of length 2")
  }

  y_c <- records$certain
  y_u <- records$uncertain

  # Calculate key values
  bigT <- length(y_c)
  t_m <- max(which(y_c > 0))
  total_c <- sum(y_c)
  total_u <- sum(y_u)
  cumsum_u <- cumsum(y_u)
  logfact_c <- sum(lfactorial(y_c))
  logfact_u <- sum(lfactorial(y_u))
  max_lag_obs <- bigT - t_m

  if (priors$g == 0) {
    max_lag_prior <- max_lag_obs
  } else {
    lag_req <- ceiling(log(1 - tol) / log(priors$g / (priors$g + 1)) - 1)
    max_lag_prior <- max(max_lag_obs, lag_req)
  }
  n_cat <- max_lag_prior + 2

  # Check if the maximum lag before swapping to tail approximation is large
  if (max_lag_prior > 1000) {
    warning(paste0(
      "maximum lag is very large (", max_lag_prior, "); JAGS may be slow"
    ))
  }

  # Specify model and parameters
  data_list <- list(
    zeros = 0L,
    init_time = init.time,
    bigT = bigT,
    t_m = t_m,
    total_c = total_c,
    total_u = total_u,
    cumsum_u = cumsum_u,
    logfact_c = logfact_c,
    logfact_u = logfact_u,
    max_lag_prior = max_lag_prior,
    n_cat = n_cat,
    g = priors$g,
    theta_max_a = priors$theta_max[1],
    theta_max_b = priors$theta_max[2],
    lambda_v_a = priors$lambda_v[1],
    lambda_v_b = priors$lambda_v[2],
    lambda_i_a = priors$lambda_i[1],
    lambda_i_b = priors$lambda_i[2]
  )

  model_string <- "
    model {
      # 1. Priors
      theta_max ~ dbeta(theta_max_a, theta_max_b)

      lambda_v ~ dgamma(lambda_v_a, lambda_v_b)
      lambda_i ~ dgamma(lambda_i_a, lambda_i_b)

      pi_e ~ dunif(0, 1)

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
      tau_E0 <- t_m + tau_L

      # 2. Likelihood
      in_window <- step(bigT - tau_E0)
      n_extant <- in_window * tau_E0 + (1 - in_window) * bigT

      total_u_e <- cumsum_u[n_extant]
      total_u_x <- total_u - total_u_e

      n_extinct <- bigT - n_extant

      mu_c <- lambda_v * pi_e
      mu_u_e <- lambda_v * (1 - pi_e) + lambda_i
      mu_u_x <- lambda_i

      loglik_c <- total_c * log(mu_c) - n_extant * mu_c - logfact_c
      loglik_u_e <- total_u_e * log(mu_u_e) - n_extant * mu_u_e
      loglik_u_x <- total_u_x * log(mu_u_x) - n_extinct * mu_u_x
      loglik <- loglik_c + loglik_u_e + loglik_u_x -logfact_u

      nll <- -loglik
      zeros ~ dpois(nll)

      tau_E <- init_time + tau_E0 - 1
    }
  "

  # Set up random seeds for JAGS that respect set.seed() etc.
  jags_seeds <- sample.int(.Machine$integer.max, n.chains, replace = FALSE)

  inits_list <- lapply(seq_len(n.chains), function(chain) {
    list(
      theta_max = runif(1, 0.01, 0.99),
      tau_cat = sample.int(n_cat, 1L),
      tau_tail = sample.int(2L * bigT + 1L, 1L) - 1L,
      lambda_v = rgamma(1, priors$lambda_v[1], priors$lambda_v[2]),
      lambda_i = rgamma(1, priors$lambda_i[1], priors$lambda_i[2]),
      pi_e = runif(1, 0.01, 0.99),
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = jags_seeds[chain]
    )
  })

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
      "tau_E", "lambda_v", "lambda_i", "pi_e", "theta_max"
    ), n.iter = n.iter, thin = n.thin)
  }))

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(samples))

  # Calculate p(extant)
  p.extant <- mean(posterior$tau_E >= test.time)

  # Calculate point estimate
  estimate <- median(posterior$tau_E)

  # Calculate credible interval bounds
  cred.int.lower <- as.numeric(quantile(posterior$tau_E, 0))
  cred.int.upper <- as.numeric(quantile(posterior$tau_E, 1 - alpha))

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    tol = tol,
    priors = priors,
    samples = samples,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(dodo_fit(output, model = "BA26B2"))
}
