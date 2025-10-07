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
                prSD = 1, alpha, PLOT = 0) {
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
  Ldens <- rep(NA, numstepsL)
  thetavals <- seq(xmax, upperlimth, length.out = numstepsth)
  thdens <- rep(NA, numstepsth)

  # Estimate Lambda

  # Increment lambda values, integrating over theta values for each
  for (i in 1:numstepsL) {
    Ldens[i] <- ifelse(Lvals[i] <= 0,
      integrate(integrand.thetasnegL, xmax, upperlimth,
        L = Lvals[i], x = x, prmean = prmean,
        prSD = prSD
      )$value,
      integrate(integrand.thetasposL, xmax, upperlimth,
        L = Lvals[i], x = x, prmean = prmean,
        prSD = prSD
      )$value
    )
  }

  # Normalize lambda pdf to unit area
  Ldens <- Ldens / sum(Ldens)

  # Calculate posterior quantities
  Lmean <- sum(Lvals * Ldens)
  Lhat <- Lmean
  Lvar <- sum(Ldens * (Lvals - Lmean)^2)

  # Estimate Theta

  # Increment theta values, integrating over lambda values for each
  for (i in 1:numstepsth) {
    thdens[i] <- (integrate(integrand.neglambdas, -Inf, 0,
      th = thetavals[i],
      x = x, prmean = prmean, prSD = prSD
    )$value +
      integrate(integrand.poslambdas, 0, Inf,
        th = thetavals[i],
        x = x, prmean = prmean, prSD = prSD
      )$value)
  }

  # Normalize theta pdf to unit area
  thdens <- thdens / sum(thdens)

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
  temp <- t(as.matrix(c(thhat, xmax, CIupper, Lhat, Lvar)))
  colnames(temp) <- c("th-hat", "xmax", "CIupper", "L-hat", "var(L)")

  return(temp)
}
