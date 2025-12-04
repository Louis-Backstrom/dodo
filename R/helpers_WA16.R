#' @title integrand.neglambdas function from Wang et al. 2016
#'
#' @description
#' Helper function. From provided code (abm38.5.r). Calculations are done on
#' log scale until the end.
#'
#' @param L vector of lambdas
#' @param th theta value (scalar)
#' @param x vector or stratigraphic positions
#' @param prmean prior mean for lambda.
#' @param prSD prior standard deviation for lambda.
#'
#' @returns a posterior
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' @noRd

integrand.neglambdas <- function(L, th, x, prmean, prSD) {
  k <- length(L)
  prior <- dnorm(L, prmean, prSD, log = TRUE) + log(1 / th)

  likelihood <- rep(NA, k)
  for (i in 1:k) {
    likelihood[i] <- sum(log((1 - L[i]) / th * 1 / (1 - x / th)^L[i]))
  }

  posterior <- prior + likelihood
  posterior <- exp(posterior)

  return(posterior)
}

integrand.neglambdas.mpfr <- function(L, th, x, prmean, prSD,
                                      precBits = 64, return.log = FALSE) {
  # Convert everything to mpfr
  L <- Rmpfr::mpfr(L, precBits)
  th <- Rmpfr::mpfr(th, precBits)
  x <- Rmpfr::mpfr(x, precBits)
  prmean <- Rmpfr::mpfr(prmean, precBits)
  prSD <- Rmpfr::mpfr(prSD, precBits)

  # Check for potential problems with integration
  # NB: original abm code "[does] not check for xmax < th; if xmax > th, an
  # error will not be generated" - this code warns instead.
  if (any(x >= th)) {
    warning("x must be smaller than th for negative lambda model")
  }

  safe_log <- function(x, precBits = 64) {
    negorzero <- Rmpfr::asNumeric(x <= 0) == 1
    if (any(negorzero)) {
      x[negorzero] <- Rmpfr::mpfr(1e-100, precBits = precBits)
    }
    return(log(x))
  }

  k <- length(L)
  log_prior <- dlognorm(L, mean = prmean, sd = prSD) + log(1 / th)

  log_likelihood <- Rmpfr::mpfr(rep(NA, k), precBits)
  for (i in 1:k) {
    term1 <- log((1 - L[i]) / th)
    term2 <- -L[i] * safe_log(1 - x / th, precBits)
    log_likelihood[i] <- sum(term1 + term2)
  }

  log_posterior <- log_prior + log_likelihood

  if (return.log == TRUE) {
    return(log_posterior)
  } else {
    posterior <- exp(log_posterior)
    return(posterior)
  }
}

#' @title integrand.poslambdas function from Wang et al. 2016
#'
#' @description
#' Helper function. From provided code (abm38.5.r). Calculations are done on
#' log scale until the end.
#'
#' @param L vector of lambdas
#' @param th theta value (scalar)
#' @param x vector or stratigraphic positions
#' @param prmean prior mean for lambda.
#' @param prSD prior standard deviation for lambda.
#'
#' @returns a posterior
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' @noRd

integrand.poslambdas <- function(L, th, x, prmean, prSD) {
  k <- length(L)
  prior <- dnorm(L, prmean, prSD, log = TRUE) + log(1 / th)

  likelihood <- rep(NA, k)
  for (i in 1:k) {
    likelihood[i] <- sum(log((1 + L[i]) / th * (x / th)^L[i]))
  }

  posterior <- prior + likelihood
  posterior <- exp(posterior)

  return(posterior)
}

integrand.poslambdas.mpfr <- function(L, th, x, prmean, prSD,
                                      precBits = 64, return.log = FALSE) {
  # Convert everything to mpfr
  L <- Rmpfr::mpfr(L, precBits)
  th <- Rmpfr::mpfr(th, precBits)
  x <- Rmpfr::mpfr(x, precBits)
  prmean <- Rmpfr::mpfr(prmean, precBits)
  prSD <- Rmpfr::mpfr(prSD, precBits)

  k <- length(L)
  log_prior <- dlognorm(L, mean = prmean, sd = prSD) + log(1 / th)

  log_likelihood <- Rmpfr::mpfr(rep(NA, k), precBits)
  for (i in 1:k) {
    term1 <- log((1 + L[i]) / th)
    term2 <- L[i] * log(x / th)
    log_likelihood[i] <- sum(term1 + term2)
  }

  log_posterior <- log_prior + log_likelihood

  if (return.log == TRUE) {
    return(log_posterior)
  } else {
    posterior <- exp(log_posterior)
    return(posterior)
  }
}

#' @title integrand.thetasnegL function from Wang et al. 2016
#'
#' @description
#' Helper function. From provided code (abm38.5.r). Calculations are done on
#' log scale until the end.
#'
#' @param th vector of thetas
#' @param L lambda value (scalar)
#' @param x vector or stratigraphic positions
#' @param prmean prior mean for lambda.
#' @param prSD prior standard deviation for lambda.
#'
#' @returns a posterior
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' @noRd

integrand.thetasnegL <- function(th, L, x, prmean, prSD) {
  k <- length(th)
  prior <- dnorm(L, prmean, prSD, log = TRUE) + log(1 / th)

  likelihood <- rep(NA, k)
  for (i in 1:k) {
    likelihood[i] <- sum(log((1 - L) / th[i] * 1 / (1 - x / th[i])^L))
  }

  posterior <- prior + likelihood
  posterior <- exp(posterior)

  return(posterior)
}

integrand.thetasnegL.mpfr <- function(th, L, x, prmean, prSD,
                                      precBits = 64, return.log = FALSE) {
  # Convert everything to mpfr
  th <- Rmpfr::mpfr(th, precBits)
  L <- Rmpfr::mpfr(L, precBits)
  x <- Rmpfr::mpfr(x, precBits)
  prmean <- Rmpfr::mpfr(prmean, precBits)
  prSD <- Rmpfr::mpfr(prSD, precBits)

  k <- length(th)
  log_prior <- dlognorm(L, mean = prmean, sd = prSD) + log(1 / th)

  log_likelihood <- Rmpfr::mpfr(rep(NA, k), precBits)
  for (i in 1:k) {
    term1 <- log((1 - L) / th[i])
    term2 <- -L * log(1 - x / th[i])
    log_likelihood[i] <- sum(term1 + term2)
  }

  log_posterior <- log_prior + log_likelihood

  if (return.log == TRUE) {
    return(log_posterior)
  } else {
    posterior <- exp(log_posterior)
    return(posterior)
  }
}

#' @title integrand.thetasposL function from Wang et al. 2016
#'
#' @description
#' Helper function. From provided code (abm38.5.r). Calculations are done on
#' log scale until the end.
#'
#' @param th vector of thetas
#' @param L lambda value (scalar)
#' @param x vector or stratigraphic positions
#' @param prmean prior mean for lambda.
#' @param prSD prior standard deviation for lambda.
#'
#' @returns a posterior
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' @noRd

integrand.thetasposL <- function(th, L, x, prmean, prSD) {
  k <- length(th)
  prior <- dnorm(L, prmean, prSD, log = TRUE) + log(1 / th)

  likelihood <- rep(NA, k)
  for (i in 1:k) {
    likelihood[i] <- sum(log((1 + L) / th[i] * (x / th[i])^L))
  }

  posterior <- prior + likelihood
  posterior <- exp(posterior)

  return(posterior)
}

integrand.thetasposL.mpfr <- function(th, L, x, prmean, prSD,
                                      precBits = 64, return.log = FALSE) {
  # Convert everything to mpfr
  th <- Rmpfr::mpfr(th, precBits)
  L <- Rmpfr::mpfr(L, precBits)
  x <- Rmpfr::mpfr(x, precBits)
  prmean <- Rmpfr::mpfr(prmean, precBits)
  prSD <- Rmpfr::mpfr(prSD, precBits)

  k <- length(th)
  log_prior <- dlognorm(L, mean = prmean, sd = prSD) + log(1 / th)

  log_likelihood <- Rmpfr::mpfr(rep(NA, k), precBits)
  for (i in 1:k) {
    term1 <- log((1 + L) / th[i])
    term2 <- L * log(x / th[i])
    log_likelihood[i] <- sum(term1 + term2)
  }

  log_posterior <- log_prior + log_likelihood

  if (return.log == TRUE) {
    return(log_posterior)
  } else {
    posterior <- exp(log_posterior)
    return(posterior)
  }
}

#' @title drefbeta function from Wang et al. 2016
#'
#' @description
#' Helper function. From provided code (abm38.5.r).
#'
#' @param x vector of quantiles
#' @param L lambda parameter
#'
#' @returns a probability density function
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' @noRd

drefbeta <- function(x, L) {
  if (L <= 0) {
    return(dbeta(x, 1, 1 - L))
  } else {
    return(dbeta(x, 1 + L, 1))
  }
}

#' @title dlognorm function
#'
#' @description
#' Helper function, for Rmpfr calculations.
#'
#' @param x value.
#' @param mean parameter of normal distribution.
#' @param sd parameter of normal distribution.
#' @param precBits number of bits of precision for Rmpfr.
#'
#' @returns a value.
#'
#' @noRd

dlognorm <- function(x, mean, sd, precBits = 64) {
  x <- Rmpfr::mpfr(x, precBits)
  mean <- Rmpfr::mpfr(mean, precBits)
  sd <- Rmpfr::mpfr(sd, precBits)

  log_term <- log(2 * pi * sd^2)
  squared_term <- (x - mean)^2 / (2 * sd^2)

  return(-0.5 * log_term - squared_term)
}

#' @title Adaptive Beta Method function from Wang et al. 2016
#'
#' @description
#' Helper function. From provided code (abm38.5.r).
#'
#' @param x vector of sighting records.
#' @param distance whether or not measurements (`x`) represent distance above a
#' base, or time (defaults to `FALSE`).
#' @param ext whether the event of interest is extinction (`TRUE`) or
#' origination (`FALSE`) time (defaults to `FALSE`).
#' @param base value to consider 0 if known; otherwise the minimum is used and
#' sample size is decreased by 1. Defaults to `init.time`
#' @param prmean prior mean for lambda. Defaults to 0.
#' @param prSD prior standard deviation for lambda. Defaults to 1.
#' @param alpha desired significance level.
#' @param PLOT whether or not to show plots (defaults to 0, 1 to show plots).
#' @param use.mpfr whether or not to use multiple-precision floating-point
#' computation via the `Rmpfr` package. Defaults to `FALSE`, but is necessary
#' for larger datasets (N > 50). Note that this is very slow, and experimental!
#' @param cores number of cores to use if using `Rmpfr` calculation (defaults
#' to `NULL`, in which case `parallel::detectCores()` is run).
#' @param pb whether to show a progress bar if using `Rmpfr` calculation.
#' Defaults to `FALSE`.
#'
#' @returns the model outputs.
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' @noRd

abm <- function(x, distance = FALSE, ext = FALSE, base = NULL, prmean = 0,
                prSD = 1, alpha, PLOT = 0, use.mpfr = FALSE, cores = NULL,
                pb = FALSE) {
  # Get confidence level
  conf <- 1 - alpha

  # Pre-process data

  # If a base is specified, check that it is valid
  if (!is.null(base)) {
    if ((distance & ext) & (base > min(x)) |
      (distance & !ext) & (base < max(x)) |
      (!distance & ext) & (base < max(x)) |
      (!distance & !ext) & (base > min(x))) {
      stop("Invalid value for base")
    }
  }

  # Convert units relative to base of section or other zero point
  xraw <- x
  if ((distance & ext) | (!distance & !ext)) {
    if (!is.null(base)) {
      x <- x - base
    }
    if (is.null(base)) {
      base <- min(x)
      x <- x - base
      x <- sort(x, decreasing = FALSE)[-1]
    }
  }
  if ((distance & !ext) | (!distance & ext)) {
    if (!is.null(base)) {
      x <- base - x
    }
    if (is.null(base)) {
      base <- max(x)
      x <- base - x
      x <- sort(x, decreasing = FALSE)[-1]
    }
  }

  # Scale data so theta is approximately 100 (for numerical stability)
  xmax <- max(x)
  n <- length(x)
  simplethhat <- (n + 1) / n * xmax
  scalefactor <- 100 / simplethhat
  x <- x * scalefactor
  xmax <- max(x)

  # Set iteration parameters
  upperlimth <- 500
  numstepsth <- 1000
  lowerlimL <- -10
  upperlimL <- 10
  numstepsL <- 40
  Lvals <- seq(lowerlimL, upperlimL, length.out = numstepsL)
  thetavals <- seq(xmax, upperlimth, length.out = numstepsth)

  if (use.mpfr == TRUE) {
    thdens <- list()
  } else {
    thdens <- rep(NA, numstepsth)
  }

  # Estimate Theta

  # Increment theta values, integrating over lambda values for each
  if (use.mpfr == TRUE) {
    # Convert constants to mpfr
    lowerlimL_mpfr <- Rmpfr::mpfr(lowerlimL, 256)
    upperlimL_mpfr <- Rmpfr::mpfr(upperlimL, 256)
    thetavals_mpfr <- Rmpfr::mpfr(thetavals, 256)
    x_mpfr <- Rmpfr::mpfr(x, 256)
    prmean_mpfr <- Rmpfr::mpfr(prmean, 256)
    prSD_mpfr <- Rmpfr::mpfr(prSD, 256)

    # Initialise parallel cluster
    num_cores <- ifelse(is.null(cores), parallel::detectCores() - 1, cores)
    cl <- parallel::makeCluster(num_cores)

    parallel::clusterEvalQ(cl, library(Rmpfr))

    parallel::clusterExport(
      cl,
      varlist = ls(envir = asNamespace("dodo")),
      envir = asNamespace("dodo")
    )

    thdens <- vector("list", numstepsth)

    # Progress bar
    if (pb == TRUE) {
      pbar <- txtProgressBar(min = 0, max = numstepsth, style = 3)
    }

    # Compute in chunks
    chunk_size <- 10
    for (start_idx in seq(1, numstepsth, by = chunk_size)) {
      end_idx <- min(start_idx + chunk_size - 1, numstepsth)
      res <- parallel::parLapplyLB(cl, start_idx:end_idx, function(i) {
        val <- Rmpfr::mpfr(
          Rmpfr::integrateR(integrand.neglambdas.mpfr,
            lowerlimL_mpfr, 0,
            th = thetavals_mpfr[i],
            x = x_mpfr, prmean = prmean_mpfr, prSD = prSD_mpfr
          )$value +
            Rmpfr::integrateR(integrand.poslambdas.mpfr,
              0, upperlimL_mpfr,
              th = thetavals_mpfr[i],
              x = x_mpfr, prmean = prmean_mpfr, prSD = prSD_mpfr
            )$value,
          precBits = 64
        )
        list(i = i, val = val)
      })
      for (r in res) {
        thdens[[r$i]] <- r$val
      }
      if (pb == TRUE) {
        setTxtProgressBar(pbar, end_idx)
      }
    }

    if (pb == TRUE) {
      close(pbar)
    }

    parallel::stopCluster(cl)
    thdens <- do.call(c, thdens)
  } else {
    for (i in 1:numstepsth) {
      thdens[i] <-
        (integrate(integrand.neglambdas, -Inf, 0,
          th = thetavals[i],
          x = x, prmean = prmean, prSD = prSD
        )$value +
          integrate(integrand.poslambdas, 0, Inf,
            th = thetavals[i],
            x = x, prmean = prmean, prSD = prSD
          )$value)
    }
  }

  # Normalize theta pdf to unit area
  thdens <- as.numeric(thdens / sum(thdens))

  # Calculate posterior quantities
  cutoff <- which.max(cumsum(thdens) >= conf)
  CIupper <- thetavals[cutoff]
  cutoff <- which.max(cumsum(thdens) >= .5)
  thmed <- thetavals[cutoff]
  thhat <- thmed

  # Collect results

  # Un-scale data back to original scale
  thhat <- thhat / scalefactor
  xmax <- xmax / scalefactor
  CIupper <- CIupper / scalefactor
  thetavals <- thetavals / scalefactor

  # Un-convert units back to original units (from units relative to base)
  if ((distance & ext) | (!distance & !ext)) {
    thhat <- thhat + base
    xmax <- xmax + base
    CIupper <- CIupper + base
    thetavals <- thetavals + base
  }
  if ((distance & !ext) | (!distance & ext)) {
    thhat <- base - thhat
    xmax <- base - xmax
    CIupper <- base - CIupper
    thetavals <- base - thetavals
  }
  temp <- t(as.matrix(c(thhat, xmax, CIupper)))
  colnames(temp) <- c("th-hat", "xmax", "CIupper")

  return(temp)
}
