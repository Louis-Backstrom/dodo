logsumexp <- function(x) {
  m <- max(x)
  return(m + log(sum(exp(x - m))))
}

logsumexp2 <- function(a, b) {
  m <- max(a, b)
  return(m + log(exp(a - m) + exp(b - m)))
}

log.integrate.stable <- function(logf, lower, upper, ..., grid.n = 100,
                                 rel.tol = 1e-8, subdivisions = 2000) {
  grid <- seq(lower, upper, length.out = grid.n)
  lg <- logf(grid, ...)
  if (!any(is.finite(lg))) {
    return(-Inf)
  }

  c <- max(lg[is.finite(lg)])

  val <- integrate(
    f = function(z) exp(logf(z, ...) - c),
    lower = lower, upper = upper,
    rel.tol = rel.tol, subdivisions = subdivisions
  )$value

  if (!is.finite(val) || val <= 0) {
    return(-Inf)
  }

  return(log(val) + c)
}

log.integrand.neglambdas <- function(L, th, x, prmean, prSD) {
  log_prior <- dnorm(L, prmean, prSD, log = TRUE) + log(1 / th)

  s_log1m <- sum(log1p(-x / th))
  n <- length(x)

  ll <- vapply(L, function(Li) {
    if (Li >= 1) {
      return(-Inf)
    }
    n * log((1 - Li) / th) - Li * s_log1m
  }, numeric(1))

  return(log_prior + ll)
}

log.integrand.poslambdas <- function(L, th, x, prmean, prSD) {
  log_prior <- dnorm(L, prmean, prSD, log = TRUE) + log(1 / th)

  s_logx <- sum(log(x / th))
  n <- length(x)

  ll <- vapply(L, function(Li) {
    n * log((1 + Li) / th) + Li * s_logx
  }, numeric(1))

  return(log_prior + ll)
}

abm <- function(x, distance = FALSE, ext = FALSE, base = NULL, prmean = 0,
                prSD = 1, alpha) {
  # Get confidence level
  conf <- 1 - alpha

  # Pre-process data

  # If a base is specified, check that it is valid
  if (!is.null(base)) {
    if ((distance & ext) & (base > min(x)) |
      (distance & !ext) & (base < max(x)) |
      (!distance & ext) & (base < max(x)) |
      (!distance & !ext) & (base > min(x))) {
      stop("invalid value for base")
    }
  }

  # Convert units relative to base of section or other zero point
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
  thetavals <- seq(xmax, upperlimth, length.out = numstepsth)

  # Estimate Theta

  # Increment theta values, integrating over lambda values for each

  # Work on log scale to avoid under/overflow
  log_thdens <- rep(-Inf, numstepsth)

  # Finite lambda bounds (adaptive, avoids +/- Inf integration)
  # 10*SD is extremely conservative for a Normal prior
  Lmax <- max(10, abs(prmean) + 10 * prSD)
  Lneg_low <- -Lmax
  Lpos_high <- Lmax

  for (i in 1:numstepsth) {
    logI_neg <- log.integrate.stable(
      logf = log.integrand.neglambdas,
      lower = Lneg_low, upper = 0,
      th = thetavals[i], x = x, prmean = prmean, prSD = prSD
    )

    logI_pos <- log.integrate.stable(
      logf = log.integrand.poslambdas,
      lower = 0, upper = Lpos_high,
      th = thetavals[i], x = x, prmean = prmean, prSD = prSD
    )

    log_thdens[i] <- logsumexp2(logI_neg, logI_pos)
  }

  logZ <- logsumexp(log_thdens)
  thdens <- exp(log_thdens - logZ)

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
