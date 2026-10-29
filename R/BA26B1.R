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
#' @param truncation.point non-negative integer giving the maximum number of
#' sampling intervals after `test.time` at which extinction can occur. Prior
#' mass beyond `test.time + truncation.point` is redistributed proportionally
#' across extinction times from `test.time` through the truncation time.
#' Probabilities before `test.time` and the overall prior probability of being
#' extant at `test.time` are unchanged.
#' @param tol tolerance for the lagged extinction hazard tail approximation. The
#' lagged hazard is exactly evaluated until it reaches `tol` of `theta_max`,
#' after which point the remaining tail is approximated using a constant hazard
#' of `theta_max`. Defaults to 0.999; if this yields slow evaluation (in
#' conjunction with a large value of `g`), consider a lower tolerance e.g. 0.99.
#' @param priors `list` with three elements: `theta_max`, `g`, and `lambda`.
#' Both `theta_max` and `lambda` are `numeric` vectors of length two, `g` is a
#' single `numeric`. The two elements in `theta_max` are the shape parameters
#' for the Beta hyperprior on \eqn{\theta_{max}}. They default to (1, 10). `g`
#' is the delay parameter for the lagged extinction hazard. It defaults to 0,
#' which recovers the constant-hazard geometric prior for extinction time. The
#' two elements in `lambda` are the shape and rate parameters for the Gamma
#' prior on \eqn{\lambda}. They default to (1, 1).
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
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' The extinction-time prior is truncated at `test.time + truncation.point`.
#' Removed prior mass is redistributed proportionally across all retained states
#' corresponding to being extant at `test.time` (i.e. extinction at or after
#' `test.time`). Consequently, the prior probability of being extant at
#' `test.time` is the same as under the untruncated prior.
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
                   truncation.point = 100, tol = 0.999,
                   priors = list(theta_max = c(1, 10), g = 0, lambda = c(1, 1)),
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

  if (!is.numeric(truncation.point) || length(truncation.point) != 1L ||
    !is.finite(truncation.point) || truncation.point < 0 ||
    truncation.point != floor(truncation.point)) {
    stop("truncation.point must be a single non-negative integer")
  }

  truncation.point <- as.integer(truncation.point)

  if (is.null(tol) || !is.numeric(tol) || length(tol) != 1L ||
    !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("tol must be a single number in (0, 1)")
  }

  # Check that records are in a valid format
  if (!is.numeric(records) || anyNA(records) || any(!is.finite(records)) ||
    any(records < 0) || any(records != floor(records))) {
    stop("records must contain only finite, non-negative integer counts")
  }

  if (length(records) < 2L) {
    stop("records must contain at least two sampling intervals")
  }

  if (!any(records > 0)) {
    stop("records must contain at least one positive count")
  }

  records <- as.integer(records)

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

  # Check that sighting-rate prior is in a valid format
  if (is.null(priors$lambda) || !is.numeric(priors$lambda) ||
    length(priors$lambda) != 2L || anyNA(priors$lambda) ||
    any(!is.finite(priors$lambda)) || any(priors$lambda <= 0)) {
    stop("priors$lambda must be a finite positive vector of length 2")
  }

  # Calculate key values
  bigT <- length(records)
  t_m <- max(which(records > 0))
  total_records <- sum(records)
  max_lag_obs <- bigT - t_m

  if (priors$g == 0) {
    max_lag_prior <- max_lag_obs
  } else {
    lag_req <- ceiling(log(1 - tol) / log(priors$g / (priors$g + 1)) - 1)
    max_lag_prior <- max(max_lag_obs, lag_req)
  }

  test_lag <- test.time - init.time - t_m + 1

  if (!is.finite(test_lag) || test_lag != floor(test_lag)) {
    stop(
      paste(
        "test.time and init.time must define an integer number",
        "of sampling intervals"
      )
    )
  }

  if (test_lag < 0) {
    stop("test.time must not precede the earliest possible extinction time")
  }

  test_lag <- as.integer(test_lag)
  max_lag <- test_lag + truncation.point
  n_cat <- max_lag + 1L

  # Check if the exact hazard section is large
  if (max_lag_prior > 1000) {
    warning(paste0(
      "maximum lag before the tail approximation is very large (",
      max_lag_prior,
      "); JAGS may be slow"
    ))
  }

  # Check if the complete truncated prior is large
  if (max_lag > 5000) {
    warning(paste0(
      "truncated extinction-time prior contains ",
      n_cat,
      " categories; JAGS may be slow"
    ))
  }

  # Specify model and parameters
  data_list <- list(
    zeros = 0L,
    init_time = init.time,
    bigT = bigT,
    t_m = t_m,
    total_records = total_records,
    max_lag_prior = max_lag_prior,
    test_lag = test_lag,
    max_lag = max_lag,
    n_cat = n_cat,
    g = priors$g,
    theta_max_a = priors$theta_max[1],
    theta_max_b = priors$theta_max[2],
    lambda_a = priors$lambda[1],
    lambda_b = priors$lambda[2]
  )

  model_string <- "
    model {
      # 1. Priors
      theta_max ~ dbeta(theta_max_a, theta_max_b)

      ## 1.1. Original extinction time prior through max_lag
      surv[1] <- 1

      for (k in 1:n_cat) {
        lag[k] <- k - 1
        theta_lag_exact[k] <- theta_max * (1 - pow(1 - 1 / (g + 1), lag[k] + 1))
        in_exact[k] <- step(max_lag_prior - lag[k] + 0.5)
        theta_lag[k] <- in_exact[k] * theta_lag_exact[k] +
          (1 - in_exact[k]) * theta_max
        p_tau[k] <- surv[k] * theta_lag[k]
        surv[k + 1] <- surv[k] * (1 - theta_lag[k])
      }

      ## 1.2. Preserve prior mass below and above test.time
      p_extant <- surv[test_lag + 1]

      p_extant_kept <- surv[test_lag + 1] - surv[max_lag + 2]
      extant_scale <- p_extant / p_extant_kept

      for (k in 1:n_cat) {
        is_extant[k] <- step(lag[k] - test_lag + 0.5)
        p_tau_adjusted[k] <- p_tau[k] *
          ((1 - is_extant[k]) + is_extant[k] * extant_scale)
      }

      tau_cat ~ dcat(p_tau_adjusted[1:n_cat])

      tau_L <- tau_cat - 1
      tau_E0 <- t_m + tau_L

      # 2. Likelihood
      in_window <- step(bigT - tau_E0)
      n_extant <- in_window * tau_E0 + (1 - in_window) * bigT

      nll <- (lambda_a + total_records) * log(lambda_b + n_extant)
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
      "tau_E", "tau_E0", "theta_max"
    ), n.iter = n.iter, thin = n.thin)
  }))

  # Reconstruct lambda from its exact conditional posterior
  samples <- coda::mcmc.list(
    lapply(samples, function(chain) {
      chain_matrix <- as.matrix(chain)

      n_extant <- pmin(chain_matrix[, "tau_E0"], bigT)

      lambda <- rgamma(
        n = nrow(chain_matrix),
        shape = priors$lambda[1] + total_records,
        rate = priors$lambda[2] + n_extant
      )

      output_matrix <- cbind(
        tau_E = chain_matrix[, "tau_E"],
        lambda = lambda,
        theta_max = chain_matrix[, "theta_max"]
      )

      chain_parameters <- coda::mcpar(chain)

      coda::mcmc(output_matrix,
        start = chain_parameters[1],
        end = chain_parameters[2], thin = chain_parameters[3]
      )
    })
  )

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
    truncation.point = truncation.point,
    max.time = test.time + truncation.point,
    tol = tol,
    priors = priors,
    samples = samples,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(dodo_fit(output, model = "BA26B1"))
}
