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
#' @param tol tolerance for the lagged extinction hazard tail approximation. The
#' lagged hazard is exactly evaluated until it reaches `tol` of `theta_max`,
#' after which point the remaining tail is approximated using a constant hazard
#' of `theta_max`. Defaults to 0.999; if this yields slow evaluation (in
#' conjunction with a large value of `g`), consider a lower tolerance e.g. 0.99.
#' @param scale.effort whether or not to scale, but not centre, each effort
#' covariate before model fitting. Defaults to `TRUE`, in which case each
#' covariate is divided by its standard deviation. Scaling can improve MCMC
#' mixing and the numerical stability of the JAGS sampler while preserving
#' meaningful zeros.
#' @param variable.selection logical. If TRUE (default), Kuo-Mallick
#' variable selection is used to account for uncertainty in which effort
#' covariates influence sighting rates. If FALSE, all effort covariates are
#' included in the model. If only a single effort covariate is supplied and
#' variable.selection = TRUE, a warning is issued and variable selection is
#' automatically disabled.
#' @param priors `list` with four elements: `theta_max`, `g`, `rate`, and
#' `inclusion`. `theta_max` and `rate` are `numeric` vectors, `inclusion` is
#' either `NULL` or a `numeric` vector, and `g` is a single `numeric`.
#' The two elements of `theta_max` are the shape parameters for the Beta
#' hyperprior on \eqn{\theta_{max}} and default to (1, 10). `g` is the delay
#' parameter for the lagged extinction hazard and defaults to 0, which recovers
#' the constant-hazard geometric prior for extinction time. `rate` should have
#' length one or the same length as the number of effort covariates, i.e.
#' `ncol(effort)`, and specifies the rate of the exponential prior on each
#' non-negative effort coefficient. A scalar value is applied to every
#' coefficient. It defaults to 10, corresponding to a prior mean and standard
#' deviation of 0.1. `inclusion` specifies the prior probability that each
#' effort covariate is included and should have length one or `ncol(effort)`. If
#' variable.selection = TRUE and inclusion is NULL, each inclusion probability
#' defaults to 1 / ncol(effort), corresponding to a prior expected model size of
#' one before conditioning on the selected covariates providing positive effort
#' in all intervals containing positive sighting counts. If
#' variable.selection = FALSE, all effort covariates are forced into the model
#' and inclusion priors are ignored.
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
#' All sighting records are assumed to be certain. The expected sighting rate
#' is a non-negative additive function of the effort covariates, with no
#' intercept. Consequently, intervals in which every effort covariate is zero
#' must also have zero sighting counts.
#'
#' If variable.selection = TRUE, uncertainty in the effort covariates included
#' in the model is represented using binary inclusion indicators. The extinction
#' time posterior is averaged over the effort models visited by the MCMC
#' sampler.
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
                   test.time = init.time + length(records) - 1, tol = 0.999,
                   scale.effort = TRUE, variable.selection = TRUE,
                   priors = list(
                     theta_max = c(1, 10), g = 0, rate = 10, inclusion = NULL
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

  if (!is.logical(scale.effort) || length(scale.effort) != 1L ||
    is.na(scale.effort)) {
    stop("scale.effort must be TRUE or FALSE")
  }

  if (!is.logical(variable.selection) || length(variable.selection) != 1L ||
    is.na(variable.selection)) {
    stop("variable.selection must be TRUE or FALSE")
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

  if (is.null(priors$g) || !is.numeric(priors$g) || length(priors$g) != 1L ||
    !is.finite(priors$g) || priors$g < 0) {
    stop("priors$g must be a single finite non-negative number")
  }

  # Check that effort data are in a valid format
  effort <- as.matrix(effort)

  if (!is.numeric(effort) || anyNA(effort) || any(!is.finite(effort))) {
    stop("effort must contain only finite numeric values")
  }

  if (ncol(effort) < 1L) {
    stop("effort must contain at least one variable")
  }

  if (nrow(effort) != length(records)) {
    stop("effort must have one row per record")
  }

  if (any(effort < 0)) {
    stop("effort variables must be non-negative")
  }

  if (is.null(colnames(effort))) {
    colnames(effort) <- paste0("effort", seq_len(ncol(effort)))
  }

  effort_original <- effort

  if (variable.selection == TRUE && ncol(effort) == 1L) {
    warning(
      paste(
        "variable.selection = TRUE but only one effort covariate was supplied;",
        "selection has been disabled and the covariate is forced into the model"
      )
    )
    variable.selection <- FALSE
  }

  # Check that effort coefficient priors are in a valid format
  if (is.null(priors$rate) || !is.numeric(priors$rate) || anyNA(priors$rate) ||
    any(!is.finite(priors$rate)) || any(priors$rate <= 0)) {
    stop("priors$rate must contain finite positive values")
  }

  if (!length(priors$rate) %in% c(1L, ncol(effort))) {
    stop("priors$rate must have length 1 or ncol(effort)")
  }

  if (length(priors$rate) == 1L) {
    rate <- rep(priors$rate, ncol(effort))
  } else {
    rate <- priors$rate
  }

  names(rate) <- colnames(effort)

  # Check that effort variable inclusion priors are in a valid format
  if (variable.selection == TRUE) {
    if (is.null(priors$inclusion)) {
      inclusion_prob <- rep(1 / ncol(effort), ncol(effort))
    } else {
      if (!is.numeric(priors$inclusion) ||
        anyNA(priors$inclusion) || any(!is.finite(priors$inclusion)) ||
        any(priors$inclusion <= 0) || any(priors$inclusion > 1)) {
        stop("priors$inclusion must contain values in (0, 1]")
      }

      if (!length(priors$inclusion) %in% c(1L, ncol(effort))) {
        stop("priors$inclusion must have length 1 or ncol(effort)")
      }

      if (length(priors$inclusion) == 1L) {
        inclusion_prob <- rep(priors$inclusion, ncol(effort))
      } else {
        inclusion_prob <- priors$inclusion
      }
    }
  } else {
    inclusion_prob <- rep(1, ncol(effort))
  }

  names(inclusion_prob) <- colnames(effort)

  # Scale effort data without centring
  effort_scale <- rep(1, ncol(effort))
  names(effort_scale) <- colnames(effort)

  if (scale.effort == TRUE) {
    effort_scale <- apply(effort, MARGIN = 2, FUN = sd)

    if (any(!is.finite(effort_scale)) || any(effort_scale <= 0)) {
      stop(
        paste(
          "Cannot scale effort covariates;",
          "all variables must have finite, non-zero variance"
        )
      )
    }

    effort <- sweep(effort, MARGIN = 2, STATS = effort_scale, FUN = "/")
  }

  # Check compatibility with the zero-effort model assumption
  zero_effort <- rowSums(effort) == 0

  if (any(zero_effort & records > 0)) {
    bad_rows <- which(zero_effort & records > 0)

    stop(
      paste0(
        "Non-zero sighting counts occur in interval(s) with zero effort: ",
        paste(bad_rows, collapse = ", ")
      )
    )
  }

  # Calculate key values
  bigT <- length(records)
  t_m <- max(which(records > 0))
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
    init_time = init.time,
    y = records,
    x = effort,
    p = ncol(effort),
    bigT = bigT,
    t_m = t_m,
    max_lag_prior = max_lag_prior,
    n_cat = n_cat,
    g = priors$g,
    theta_max_a = priors$theta_max[1],
    theta_max_b = priors$theta_max[2],
    rate = rate,
    inclusion_prob = inclusion_prob
  )

  model_string <- "
    model {
      # 1. Priors
      theta_max ~ dbeta(theta_max_a, theta_max_b)

      for (m in 1:p) {
        included[m] ~ dbern(inclusion_prob[m])
        beta_raw[m] ~ dexp(rate[m])
        beta[m] <- included[m] * beta_raw[m]
      }

      n_included <- sum(included[1:p])

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
      for (t in 1:bigT) {
        lambda[t] <- inprod(beta[1:p], x[t, 1:p])

        extant[t] <- step(tau_E0 - t + 0.5)

        mu[t] <- extant[t] * lambda[t]

        y[t] ~ dpois(mu[t])
      }

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
      included = rep(1L, ncol(effort)),
      beta_raw = pmax(rexp(ncol(effort), rate = rate), 1e-6),
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
    variables <- c("tau_E", "theta_max", "beta")

    if (variable.selection == TRUE) {
      variables <- c(variables, "beta_raw", "included", "n_included")
    }

    samples <- rjags::coda.samples(
      jags_model,
      variable.names = variables,
      n.iter = n.iter,
      thin = n.thin
    )
  }))

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(samples))

  # Calculate posterior inclusion probabilities
  if (variable.selection == TRUE) {
    included_columns <- grep(
      "^included($|\\[)",
      colnames(posterior),
      value = TRUE
    )

    if (length(included_columns) != ncol(effort)) {
      stop("could not identify all effort variable inclusion samples")
    }

    posterior_inclusion <- colMeans(
      posterior[, included_columns, drop = FALSE]
    )

    names(posterior_inclusion) <- colnames(effort_original)

    included_samples <- posterior[, included_columns, drop = FALSE]

    model_id <- apply(included_samples, MARGIN = 1, FUN = paste0, collapse = "")

    posterior_models <- sort(prop.table(table(model_id)), decreasing = TRUE)

    model_variables <- vapply(
      strsplit(names(posterior_models), split = ""),
      function(x) {
        selected <- as.logical(as.integer(x))
        paste(colnames(effort_original)[selected], collapse = " + ")
      },
      character(1)
    )

    names(posterior_models) <- model_variables
  } else {
    posterior_inclusion <- setNames(
      rep(1, ncol(effort_original)), colnames(effort_original)
    )

    posterior_models <- structure(
      1,
      names = paste(colnames(effort_original), collapse = " + ")
    )
  }

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
    effort = effort_original,
    effort.model = effort,
    effort.scale = effort_scale,
    inclusion.prob = inclusion_prob,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    tol = tol,
    scale.effort = scale.effort,
    priors = priors,
    samples = samples,
    posterior.inclusion = posterior_inclusion,
    posterior.models = posterior_models,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(dodo_fit(output, model = "BA26B3"))
}
