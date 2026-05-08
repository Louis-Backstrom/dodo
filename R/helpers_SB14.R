sb14.extended.model <- function(DATA, inputs, modelnumber, gamma,
                                increment = 0.01, increment2 = 0.01,
                                SO12 = FALSE) {
  data <- DATA$data
  indmc <- (data[, 2] == 1)
  tMin <- min(data[, 1])
  tCert <- max(data[indmc, 1])
  if (sum(indmc) == 0) {
    tCert <- tMin
  }

  DATA$tCert <- tCert
  DATA$tMin <- tMin
  n <- nrow(data)
  DATA$N <- n
  DATA$logfact <- c(0, cumsum(log(seq_len(n))))

  tmp <- sub.estimate.rate(DATA, increment, inputs)
  rate <- tmp$rate
  rateV <- tmp$rateV

  if (tCert == inputs$T) {
    BayesFactor <- 0
    warning(paste(
      "last certain sighting is at the end of the sighting period:",
      "setting Bayes Factor to 0"
    ))
    return(list(DATA = DATA, Results = BayesFactor, rate = rate))
  }

  tL <- (tCert - tMin) / (inputs$T - tMin)
  tL <- round((1 / increment2) * tL) * increment2

  if (tL + increment > 1 - (increment / 2)) {
    stop("tL is too close to 1: try making increment smaller")
  }

  tLV <- seq(from = (tL + increment), to = 1 - (increment / 2), by = increment2)
  dVec <- (data[, 1] - tMin) / (inputs$T - tMin)

  if ((modelnumber != 1) & (modelnumber != 2)) {
    stop("`modelnumber` not defined")
  } else {
    p_t_tau_E <- likelihood.all(
      tLV = tLV, indmc = indmc, tL = tL,
      increment = increment, n = n, dVec = dVec, modelnumber = modelnumber,
      SO12 = SO12, logfact = DATA$logfact
    )
  }

  if (is.null(p_t_tau_E)) {
    Results <- -Inf
  }

  prior <- as.list(NULL)
  prior$Exp <- (gamma * exp(-gamma * tLV) / (1 - exp(-gamma)))
  prior$Unif <- rep(1, length(tLV))
  prior$Tri <- (-2 * tLV) + 2

  k <- match(inputs$posteriorPrior, names(prior))
  {
    if (is.na(k)) {
      stop("no such prior")
    } else {
      p_tau_E <- prior[[k]]
    }
  }

  f <- p_t_tau_E * p_tau_E
  p_t_E <- sum(f) * increment2
  p_tauE_t <- f / p_t_E
  p_t_notE <- p_t_tau_E[length(p_t_tau_E)]

  if (p_t_notE == 0) {
    BayesFactor <- if (p_t_E == 0) NA_real_ else Inf
  } else {
    BayesFactor <- as.numeric(p_t_E / p_t_notE)
  }

  BayesFactor <- as.numeric(p_t_E / p_t_notE)

  return(list(DATA = DATA, Results = BayesFactor, rate = rate))
}

sub.estimate.rate <- function(DATA, increment, inputs) {
  icc <- 0
  rateV <- (DATA$tMin):(inputs$T)
  rate <- array(dim = c(length(rateV), 3))

  for (icc in 1:length(rateV)) {
    TE <- rateV[icc]
    indm <- (DATA$data[, 1] <= TE)
    m <- sum(indm)
    yr <- TE
    numBefore <- m
    timeBefore <- yr - DATA$tMin
    numAfter <- DATA$N - numBefore
    timeAfter <- inputs$T - yr
    rate[icc, 1] <- numBefore / timeBefore
    rate[icc, 2] <- numAfter / timeAfter
    rate[icc, 3] <- yr
  }
  out <- list(rate = rate, rateV = rateV)
  return(out)
}

likelihood.all <- function(tLV, indmc, tL, increment, n, dVec,
                           modelnumber, SO12 = FALSE, logfact) {
  log_integrand <- rep(-Inf, length(tLV))

  for (iy in seq_along(tLV)) {
    TE <- tLV[iy]
    indm <- (dVec <= TE)
    m <- sum(indm)
    mc <- sum(indmc * indm)

    if (TE > tL) {
      log_integrand[iy] <- likelihood.integrated(
        TE = TE, n = n, m = m, mc = mc,
        modelnumber = modelnumber,
        increment = increment, tL = tL,
        SO12 = SO12, logfact = logfact
      )
    }
  }

  # numerically safe exponentiation by shifting
  finite <- is.finite(log_integrand)
  if (!any(finite)) {
    return(rep(0, length(tLV)))
  }

  M <- max(log_integrand[finite])
  integrand <- numeric(length(tLV))
  integrand[finite] <- exp(log_integrand[finite] - M)
  integrand[!finite] <- 0

  return(integrand)
}

likelihood.integrated <- function(TE, n, m, mc, modelnumber, increment,
                                  tL, SO12 = FALSE, logfact) {
  omega <- seq(
    from = (increment / 2),
    to = (1 - (increment / 2)),
    by = increment
  )

  log_vals <- numeric(length(omega))

  for (i in seq_along(omega)) {
    log_vals[i] <- likelihood.each(
      omega = omega[i], TE = TE, n = n,
      m = m, mc = mc, modelnumber = modelnumber,
      tL = tL, SO12 = SO12, logfact = logfact
    )
  }

  maxlv <- max(log_vals)
  log_out <- maxlv + log(sum(exp(log_vals - maxlv))) + log(increment)
  return(log_out)
}

likelihood.each <- function(omega, TE, n, m, mc, modelnumber, tL,
                            SO12 = FALSE, logfact) {
  z <- (1 - omega) / omega

  if (modelnumber == 2) {
    part0 <- logfact[mc + 1] +
      logfact[n - mc + 1] +
      mc * log(1 / TE) +
      mc * log(TE + z)
  } else {
    part0 <- 0
  }

  j <- mc:m

  termj1 <- logfact[m - mc + 1] -
    logfact[j - mc + 1] -
    logfact[m - j + 1]

  termj2 <- (n - j) * log(z)

  if (SO12 == TRUE) {
    termj3 <- -n * log(TE + z * (1 - tL))
  } else {
    termj3 <- -n * log(TE + z)
  }

  log_terms <- part0 + termj1 + termj2 + termj3

  maxlt <- max(log_terms)
  log_out <- maxlt + log(sum(exp(log_terms - maxlt)))

  return(log_out)
}
