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
#' variable selection is to account for uncertainty in which effort covariates
#' influence sighting rates of both valid and invalid sightings. If FALSE, all
#' effort covariates are included in both submodels. If only a single effort
#' covariate is supplied and variable.selection = TRUE, a warning is issued and
#' variable selection is automatically disabled.
#' @param priors `list` with six elements: `theta_max`, `g`, `rate_v`,
#' `rate_i`, `inclusion_v`, and `inclusion_i`. `theta_max`, `rate_v`, and
#' `rate_i` are `numeric` vectors; `inclusion_v` and `inclusion_i` are either
#' `NULL` or `numeric` vectors; and `g` is a single `numeric`. The two elements
#' of `theta_max` are the shape parameters for the Beta hyperprior on
#' \eqn{\theta_{max}} and default to (1, 10). `g` is the delay parameter for
#' the lagged extinction hazard and defaults to 0, which recovers the
#' constant-hazard geometric prior. `rate_v` and `rate_i` specify the rates of
#' the exponential priors on the non-negative valid- and invalid-sighting
#' effort coefficients, respectively. Each must have length one or
#' `ncol(effort)` and defaults to 10, corresponding to a prior mean and
#' standard deviation of 0.1. `inclusion_v` and `inclusion_i` specify the prior
#' probabilities that each effort covariate is included in the valid- and
#' invalid-sighting submodels. If variable.selection = TRUE and inclusion_v or
#' inclusion_i is NULL, each inclusion probability defaults to 1 / ncol(effort).
#' If variable.selection = FALSE, all effort covariates are forced into both
#' submodels and inclusion priors are ignored.
#' @param n.chains number of MCMC chains to run. Defaults to 4.
#' @param n.iter number of iterations in each chain. Defaults to 110,000.
#' @param n.burnin number of iterations to discard as burn-in. Defaults to
#' 10,000.
#' @param n.thin thinning rate. Defaults to 10.
#'
#' @returns a `dodo_fit` object containing the original sighting and effort
#' data, the scaled effort matrix used for fitting, effort scaling factors,
#' posterior samples, posterior inclusion and model probabilities for both
#' sighting submodels, p(extant), the point estimate, and the credible interval.
#'
#' @note
#' This model incorporates both sighting uncertainty and variable survey effort.
#' The expected valid- and invalid-sighting rates are non-negative additive
#' functions of the effort covariates, with no intercepts. Consequently,
#' intervals in which every effort covariate is zero must also have zero
#' certain and uncertain sighting counts.
#'
#' If variable.selection = TRUE, uncertainty in the effort covariates included
#' in the valid- and invalid-sighting submodels is represented using separate
#' binary inclusion indicators. Extinction time summaries are therefore averaged
#' over the effort models visited by the MCMC sampler.
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
                   test.time = init.time + nrow(records) - 1, tol = 0.999,
                   scale.effort = TRUE, variable.selection = TRUE,
                   priors = list(
                     theta_max = c(1, 10), g = 0,
                     rate_v = 10, rate_i = 10,
                     inclusion_v = NULL, inclusion_i = NULL
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

  if (!is.logical(variable.selection) ||
    length(variable.selection) != 1L ||
    is.na(variable.selection)) {
    stop("variable.selection must be TRUE or FALSE")
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

  # Check that effort data are in a valid format
  effort <- as.matrix(effort)

  if (!is.numeric(effort) || anyNA(effort) ||
    any(!is.finite(effort))) {
    stop("effort must contain only finite numeric values")
  }

  if (ncol(effort) < 1L) {
    stop("effort must contain at least one variable")
  }

  if (nrow(effort) != nrow(records)) {
    stop("effort must have one row per record")
  }

  if (any(effort < 0)) {
    stop("effort variables must be non-negative")
  }

  if (is.null(colnames(effort))) {
    colnames(effort) <- paste0("effort", seq_len(ncol(effort)))
  }

  effort_original <- effort

  if (variable.selection && ncol(effort) == 1L) {
    warning(
      paste(
        "variable.selection = TRUE but only one effort covariate was supplied;",
        "selection has been disabled and the covariate will be forced into",
        "both sighting submodels"
      )
    )

    variable.selection <- FALSE
  }

  # Check that effort coefficient priors are in a valid format
  if (is.null(priors$rate_v) || !is.numeric(priors$rate_v) ||
    anyNA(priors$rate_v) || any(!is.finite(priors$rate_v)) ||
    any(priors$rate_v <= 0)) {
    stop("priors$rate_v must contain finite positive values")
  }

  if (is.null(priors$rate_i) || !is.numeric(priors$rate_i) ||
    anyNA(priors$rate_i) || any(!is.finite(priors$rate_i)) ||
    any(priors$rate_i <= 0)) {
    stop("priors$rate_i must contain finite positive values")
  }

  if (!length(priors$rate_v) %in% c(1L, ncol(effort))) {
    stop("priors$rate_v must have length 1 or ncol(effort)")
  }

  if (!length(priors$rate_i) %in% c(1L, ncol(effort))) {
    stop("priors$rate_i must have length 1 or ncol(effort)")
  }

  if (length(priors$rate_v) == 1L) {
    rate_v <- rep(priors$rate_v, ncol(effort))
  } else {
    rate_v <- priors$rate_v
  }

  if (length(priors$rate_i) == 1L) {
    rate_i <- rep(priors$rate_i, ncol(effort))
  } else {
    rate_i <- priors$rate_i
  }

  names(rate_v) <- colnames(effort)
  names(rate_i) <- colnames(effort)

  # Check that effort variable inclusion priors are in a valid format
  if (variable.selection == TRUE) {
    if (is.null(priors$inclusion_v)) {
      inclusion_prob_v <- rep(1 / ncol(effort), ncol(effort))
    } else {
      if (!is.numeric(priors$inclusion_v) ||
        anyNA(priors$inclusion_v) ||
        any(!is.finite(priors$inclusion_v)) ||
        any(priors$inclusion_v <= 0) ||
        any(priors$inclusion_v > 1)) {
        stop("priors$inclusion_v must contain values in (0, 1]")
      }

      if (!length(priors$inclusion_v) %in% c(1L, ncol(effort))) {
        stop("priors$inclusion_v must have length 1 or ncol(effort)")
      }

      if (length(priors$inclusion_v) == 1L) {
        inclusion_prob_v <- rep(priors$inclusion_v, ncol(effort))
      } else {
        inclusion_prob_v <- priors$inclusion_v
      }
    }
  } else {
    inclusion_prob_v <- rep(1, ncol(effort))
  }

  if (variable.selection == TRUE) {
    if (is.null(priors$inclusion_i)) {
      inclusion_prob_i <- rep(1 / ncol(effort), ncol(effort))
    } else {
      if (!is.numeric(priors$inclusion_i) ||
        anyNA(priors$inclusion_i) ||
        any(!is.finite(priors$inclusion_i)) ||
        any(priors$inclusion_i <= 0) ||
        any(priors$inclusion_i > 1)) {
        stop("priors$inclusion_i must contain values in (0, 1]")
      }

      if (!length(priors$inclusion_i) %in% c(1L, ncol(effort))) {
        stop("priors$inclusion_i must have length 1 or ncol(effort)")
      }

      if (length(priors$inclusion_i) == 1L) {
        inclusion_prob_i <- rep(priors$inclusion_i, ncol(effort))
      } else {
        inclusion_prob_i <- priors$inclusion_i
      }
    }
  } else {
    inclusion_prob_i <- rep(1, ncol(effort))
  }

  names(inclusion_prob_v) <- colnames(effort)
  names(inclusion_prob_i) <- colnames(effort)

  # Scale effort data without centring
  effort_scale <- rep(1, ncol(effort))
  names(effort_scale) <- colnames(effort)

  if (scale.effort) {
    effort_scale <- apply(effort, MARGIN = 2, FUN = stats::sd)

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
  positive_sightings <- records$certain > 0 | records$uncertain > 0

  if (any(zero_effort & positive_sightings)) {
    bad_rows <- which(zero_effort & positive_sightings)

    stop(
      paste0(
        "Non-zero sighting counts occur in interval(s) with zero effort: ",
        paste(bad_rows, collapse = ", ")
      )
    )
  }

  # Identify effort availability and positive-sighting intervals
  effort_positive <- 1L * (effort > 0)
  certain_positive <- as.integer(records$certain > 0)
  uncertain_positive <- as.integer(records$uncertain > 0)

  y_c <- records$certain
  y_u <- records$uncertain

  # Calculate key values
  bigT <- nrow(records)
  t_m <- max(which(y_c > 0))
  max_lag_obs <- bigT - t_m
  haz_tol <- tol

  if (priors$g == 0) {
    max_lag_prior <- max_lag_obs
  } else {
    lag_req <- ceiling(
      log(1 - haz_tol) / log(priors$g / (priors$g + 1)) - 1
    )

    max_lag_prior <- max(max_lag_obs, lag_req)
  }

  n_cat <- max_lag_prior + 2
  y_c_logfact <- lfactorial(y_c)
  y_u_logfact <- lfactorial(y_u)

  # Check if the maximum lag before swapping to tail approximation is large
  if (max_lag_prior > 1000) {
    warning(paste0(
      "maximum lag is very large (", max_lag_prior, "); JAGS may be slow"
    ))
  }

  # Specify model and parameters
  data_list <- list(
    init_time = init.time,
    y_c = y_c,
    y_u = y_u,
    x = effort,
    p = ncol(effort),
    bigT = bigT,
    t_m = t_m,
    max_lag_prior = max_lag_prior,
    n_cat = n_cat,
    y_c_logfact = y_c_logfact,
    y_u_logfact = y_u_logfact,
    zeros = 0L,
    valid_model = 1L,
    effort_positive = effort_positive,
    certain_positive = certain_positive,
    uncertain_positive = uncertain_positive,
    g = priors$g,
    theta_max_a = priors$theta_max[1],
    theta_max_b = priors$theta_max[2],
    rate_v = rate_v,
    rate_i = rate_i,
    inclusion_prob_v = inclusion_prob_v,
    inclusion_prob_i = inclusion_prob_i
  )

  model_string <- "
    model {
      # 1. Priors
      theta_max ~ dbeta(theta_max_a, theta_max_b)

      pi_e ~ dunif(0, 1)

      for (m in 1:p) {
        included_v[m] ~ dbern(inclusion_prob_v[m])
        included_i[m] ~ dbern(inclusion_prob_i[m])

        beta_raw[m] ~ dexp(rate_v[m])
        gamma_raw[m] ~ dexp(rate_i[m])

        beta[m] <- included_v[m] * beta_raw[m]
        gamma[m] <- included_i[m] * gamma_raw[m]
      }

      n_included_v <- sum(included_v[1:p])
      n_included_i <- sum(included_i[1:p])

      ## 1.1. Time-varying hazard prior for tau_L up to max_lag_prior
      surv[1] <- 1

      for (k in 1:(max_lag_prior + 1)) {
        lag[k] <- k - 1
        theta_lag[k] <- theta_max *
          (1 - pow(1 - 1 / (g + 1), lag[k] + 1))
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

      # Ensure selected covariates are compatible with observed sightings
      for (t in 1:bigT) {
        selected_effort_v[t] <- inprod(
          included_v[1:p],
          effort_positive[t, 1:p]
        )

        selected_effort_i[t] <- inprod(
          included_i[1:p],
          effort_positive[t, 1:p]
        )

        effort_v_covered[t] <- step(selected_effort_v[t] - 0.5)
        effort_i_covered[t] <- step(selected_effort_i[t] - 0.5)

        effort_any_covered[t] <- step(
          selected_effort_v[t] + selected_effort_i[t] - 0.5
        )

        # Extinction at tau_E0 means periods through tau_E0 are extant
        post_extinction[t] <- step(t - tau_E0 - 0.5)

        # Certain sightings require valid-sighting effort
        valid_certain[t] <- 1 - certain_positive[t] *
          (1 - effort_v_covered[t])

        # Before extinction, uncertain sightings may be valid or invalid.
        # After extinction, uncertain sightings must be invalid.
        valid_uncertain_pre[t] <- 1 - uncertain_positive[t] *
          (1 - effort_any_covered[t])

        valid_uncertain_post[t] <- 1 - uncertain_positive[t] *
          (1 - effort_i_covered[t])

        valid_uncertain[t] <- (1 - post_extinction[t]) *
          valid_uncertain_pre[t] +
          post_extinction[t] * valid_uncertain_post[t]

        valid_period[t] <- valid_certain[t] * valid_uncertain[t]
      }

      all_periods_valid <- prod(valid_period[1:bigT])
      valid_model ~ dbern(all_periods_valid)

      # 2. Likelihood
      for (t in 1:bigT) {
        lambda_v[t] <- inprod(beta[1:p], x[t, 1:p])
        lambda_i[t] <- inprod(gamma[1:p], x[t, 1:p])

        mu_c_ext[t] <- lambda_v[t] * pi_e
        mu_u_ext[t] <- lambda_v[t] * (1 - pi_e) + lambda_i[t]

        mu_c_safe[t] <- mu_c_ext[t] + equals(mu_c_ext[t], 0)
        mu_u_safe[t] <- mu_u_ext[t] + equals(mu_u_ext[t], 0)
        lambda_i_safe[t] <- lambda_i[t] + equals(lambda_i[t], 0)

        loglik_c_ext[t] <- -mu_c_ext[t] +
          y_c[t] * log(mu_c_safe[t]) -
          y_c_logfact[t]

        loglik_u_ext[t] <- -mu_u_ext[t] +
          y_u[t] * log(mu_u_safe[t]) -
          y_u_logfact[t]

        loglik_ext[t] <- loglik_c_ext[t] + loglik_u_ext[t]

        loglik_u_post[t] <- -lambda_i[t] +
          y_u[t] * log(lambda_i_safe[t]) -
          y_u_logfact[t]
      }

      cum_ext[1] <- loglik_ext[1]
      cum_u_post[1] <- loglik_u_post[1]

      for (t in 2:bigT) {
        cum_ext[t] <- cum_ext[t - 1] + loglik_ext[t]
        cum_u_post[t] <- cum_u_post[t - 1] + loglik_u_post[t]
      }

      for (t in 1:bigT) {
        post_u_after[t] <- cum_u_post[bigT] - cum_u_post[t]
        loglik[t] <- cum_ext[t] + post_u_after[t]
      }

      loglik[bigT + 1] <- cum_ext[bigT]

      idx <- step(bigT - tau_E0) * tau_E0 +
        step(tau_E0 - bigT - 1) * (bigT + 1)

      selected_loglik <- loglik[idx]

      phi <- -selected_loglik
      zeros ~ dpois(phi)

      tau_E <- init_time + tau_E0 - 1
    }
  "

  inits_list <- function() {
    list(
      theta_max = runif(1, 0.01, 0.99),
      tau_cat = sample(1:n_cat, 1),
      tau_tail = sample(0:(2 * bigT), 1),
      pi_e = runif(1, 0.01, 0.99),
      included_v = rep(1L, ncol(effort)),
      included_i = rep(1L, ncol(effort)),
      beta_raw = pmax(
        rexp(ncol(effort), rate = rate_v),
        1e-6
      ),
      gamma_raw = pmax(
        rexp(ncol(effort), rate = rate_i),
        1e-6
      )
    )
  }

  model_file <- tempfile(fileext = ".txt")
  writeLines(model_string, con = model_file)
  on.exit(unlink(model_file), add = TRUE)

  # Run MCMC sampling
  invisible(capture.output({
    jags_model <- rjags::jags.model(
      file = model_file,
      data = data_list,
      inits = inits_list,
      n.chains = n.chains,
      n.adapt = n.burnin
    )

    update(jags_model, n.iter = n.burnin)

    variables <- c(
      "tau_E", "pi_e", "theta_max",
      "beta", "gamma"
    )

    if (variable.selection) {
      variables <- c(
        variables,
        "beta_raw", "gamma_raw",
        "included_v",
        "included_i",
        "n_included_v",
        "n_included_i"
      )
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

  # Identify posterior inclusion samples
  if (variable.selection) {
    included_columns_v <- grep(
      "^included_v($|\\[)",
      colnames(posterior),
      value = TRUE
    )

    included_columns_i <- grep(
      "^included_i($|\\[)",
      colnames(posterior),
      value = TRUE
    )

    if (length(included_columns_v) != ncol(effort)) {
      stop("could not identify all valid-sighting inclusion samples")
    }

    if (length(included_columns_i) != ncol(effort)) {
      stop("could not identify all invalid-sighting inclusion samples")
    }

    # Calculate posterior inclusion probabilities
    posterior_inclusion_v <- colMeans(
      posterior[, included_columns_v, drop = FALSE]
    )

    posterior_inclusion_i <- colMeans(
      posterior[, included_columns_i, drop = FALSE]
    )

    names(posterior_inclusion_v) <- colnames(effort_original)
    names(posterior_inclusion_i) <- colnames(effort_original)

    # Calculate posterior model probabilities
    included_samples_v <- posterior[
      ,
      included_columns_v,
      drop = FALSE
    ]

    included_samples_i <- posterior[
      ,
      included_columns_i,
      drop = FALSE
    ]

    model_id_v <- apply(
      included_samples_v,
      MARGIN = 1,
      FUN = paste0,
      collapse = ""
    )

    model_id_i <- apply(
      included_samples_i,
      MARGIN = 1,
      FUN = paste0,
      collapse = ""
    )

    model_id <- paste0(
      "valid:", model_id_v,
      "|invalid:", model_id_i
    )

    posterior_models <- sort(
      prop.table(table(model_id)),
      decreasing = TRUE
    )

    model_variables <- vapply(
      names(posterior_models),
      function(x) {
        parts <- strsplit(
          x,
          split = "\\|",
          fixed = FALSE
        )[[1]]

        valid_id <- sub("valid:", "", parts[1], fixed = TRUE)
        invalid_id <- sub("invalid:", "", parts[2], fixed = TRUE)

        selected_v <- strsplit(valid_id, split = "")[[1]] == "1"
        selected_i <- strsplit(invalid_id, split = "")[[1]] == "1"

        valid_names <- paste(
          colnames(effort_original)[selected_v],
          collapse = " + "
        )

        invalid_names <- paste(
          colnames(effort_original)[selected_i],
          collapse = " + "
        )

        paste0(
          "valid: ", valid_names,
          " | invalid: ", invalid_names
        )
      },
      character(1)
    )

    names(posterior_models) <- model_variables
  } else {
    posterior_inclusion_v <- setNames(
      rep(1, ncol(effort_original)),
      colnames(effort_original)
    )

    posterior_inclusion_i <- setNames(
      rep(1, ncol(effort_original)),
      colnames(effort_original)
    )

    posterior_models <- structure(
      1,
      names = paste0(
        "valid: ",
        paste(colnames(effort_original), collapse = " + "),
        " | invalid: ",
        paste(colnames(effort_original), collapse = " + ")
      )
    )
  }

  # Calculate p(extant)
  p.extant <- mean(posterior$tau_E >= test.time)

  # Calculate point estimate
  estimate <- median(posterior$tau_E)

  # Calculate credible interval bounds
  cred.int.lower <- as.numeric(quantile(posterior$tau_E, 0))
  cred.int.upper <- as.numeric(quantile(
    posterior$tau_E,
    1 - alpha
  ))

  # Output
  output <- list(
    records = records,
    effort = effort_original,
    effort.model = effort,
    effort.scale = effort_scale,
    inclusion.prob.v = inclusion_prob_v,
    inclusion.prob.i = inclusion_prob_i,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    tol = tol,
    scale.effort = scale.effort,
    priors = priors,
    samples = samples,
    posterior.inclusion.v = posterior_inclusion_v,
    posterior.inclusion.i = posterior_inclusion_i,
    posterior.models = posterior_models,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(dodo_fit(output, model = "BA26B4"))
}
