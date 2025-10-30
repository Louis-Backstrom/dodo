#' @title sb14.extended.model function from Solow & Beet 2014
#'
#' @description
#' Helper function. Adapted from code (woodpecker-functions-version2.r)
#' provided by Adam Butler.
#'
#' @param DATA a `list` containing all the model data.
#' @param inputs a `list` of parameters.
#' @param modelnumber which model (1 or 2) to run.
#' @param gamma a `numeric` specifying the rate parameter for the exponential
#' prior, if it is used.
#' @param increment a `numeric`.
#' @param increment2 a `numeric`.
#'
#' @returns a `list` with the model outputs.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @noRd

sb14.extended.model <- function(DATA, inputs, modelnumber, gamma,
                                increment = 0.01, increment2 = 0.01) {
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

  tL <- (tCert - tMin) / (inputs$T - tMin)
  tL <- round((1 / increment2) * tL) * increment2
  tLV <- seq(from = (tL + increment), to = 1 - (increment / 2), by = increment2)
  dVec <- (data[, 1] - tMin) / (inputs$T - tMin)

  tmp <- sub.estimate.rate(DATA, increment, inputs)
  rate <- tmp$rate
  rateV <- tmp$rateV

  if ((modelnumber != 1) & (modelnumber != 2)) {
    stop("modelnumber not defined")
  } else {
    p_t_tau_E <- likelihood.all(
      tLV = tLV, indmc = indmc, tL = tL,
      increment = increment, n = n, dVec = dVec, modelnumber = modelnumber
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
      stop("No such prior!")
    } else {
      p_tau_E <- prior[[k]]
    }
  }

  f <- p_t_tau_E * p_tau_E
  p_t_E <- sum(f) * increment2
  p_tauE_t <- f / p_t_E
  p_t_notE <- p_t_tau_E[length(p_t_tau_E)]

  BayesFactor <- as.numeric(p_t_E / p_t_notE)

  return(list(DATA = DATA, Results = BayesFactor, rate = rate))
}

#' @title sub.estimate.rate function from Solow & Beet 2014
#'
#' @description
#' Helper function. Adapted from code (woodpecker-functions-version2.r)
#' provided by Adam Butler.
#'
#' @param DATA a `list` containing all the model data.
#' @param increment a `numeric`.
#' @param inputs a `list` of parameters.
#'
#' @returns a `list` with the two rates included as elements.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @noRd

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

#' @title likelihood.all function from Solow & Beet 2014
#'
#' @description
#' Helper function. Adapted from code (woodpecker-functions-version2.r)
#' provided by Adam Butler.
#'
#' @param tLV a `numeric` vector.
#' @param indmc a `logical` vector.
#' @param tL a `numeric`.
#' @param increment a `numeric`.
#' @param n an `integer`.
#' @param dVec a `numeric` vector.
#' @param modelnumber which model (1 or 2) to run.
#'
#' @returns A `numeric` representing the likelihood integrand evaluated at each
#' possible extinction time.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @noRd

likelihood.all <- function(tLV, indmc, tL, increment, n, dVec, modelnumber) {
  integrand <- list()

  for (iy in 1:length(tLV)) {
    TE <- tLV[iy]
    indm <- (dVec <= TE)
    m <- sum(indm)
    mc <- sum(indmc * indm)

    if (TE > tL) {
      integrand[[iy]] <- likelihood.integrated(
        TE = TE, n = n, m = m, mc = mc,
        modelnumber = modelnumber, increment = increment
      )
    }
  }

  return(do.call(c, integrand))
}

#' @title likelihood.integrated function from Solow & Beet 2014
#'
#' @description
#' Helper function. Adapted from code (woodpecker-functions-version2.r)
#' provided by Adam Butler.
#'
#' @param TE a `numeric`.
#' @param n an integer.
#' @param m an integer.
#' @param mc an integer.
#' @param modelnumber which model (1 or 2) to run.
#' @param increment a `numeric`.
#'
#' @returns a value representing the integrated likelihood.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @noRd

likelihood.integrated <- function(TE, n, m, mc, modelnumber, increment) {
  omega <- seq(
    from = (increment / 2),
    to = (1 - (increment / 2)),
    by = increment
  )

  out <- 0

  for (i in 1:length(omega)) {
    out <- out + increment * likelihood.each(
      omega = omega[i], TE = TE, n = n,
      m = m, mc = mc, modelnumber = modelnumber
    )
  }

  return(out)
}

#' @title likelihood.each function from Solow & Beet 2014
#'
#' @description
#' Helper function. Adapted from code (woodpecker-functions-version2.r)
#' provided by Adam Butler.
#'
#' @param omega a `numeric` between 0 and 1.
#' @param TE a `numeric`.
#' @param n an integer.
#' @param m an integer.
#' @param mc an integer.
#' @param modelnumber which model (1 or 2) to run.
#'
#' @returns A `numeric` representing the computed likelihood.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @noRd

likelihood.each <- function(omega, TE, n, m, mc, modelnumber) {
  out <- 0
  z <- (1 - omega) / omega

  term01 <- logfact(mc - 1)
  term02 <- logfact(n - mc - 1)
  term03 <- mc * log(1 / TE)
  term04 <- mc * log(TE + z)
  part0 <- (modelnumber == 2) * (term01 + term02 + term03 + term04)

  for (j in mc:m) {
    termj1 <- logcomb(j - mc, m - mc)
    termj2 <- (n - j) * log(z)
    termj3 <- -n * log(TE + z)

    partj <- termj1 + termj2 + termj3
    new <- exp(
      Rmpfr::mpfr(part0, precBits = 64) +
        Rmpfr::mpfr(partj, precBits = 64)
    )
    out <- out + new
  }

  return(out)
}

#' @title logcomb function from Solow & Beet 2014
#'
#' @description
#' Helper function. Adapted from code (woodpecker-functions-version2.r)
#' provided by Adam Butler.
#'
#' @param a a non-negative integer.
#' @param b a non-negative integer.
#'
#' @returns a `numeric`.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @noRd

logcomb <- function(a, b) {
  if (a > b) {
    stop("Error! a > b")
  } else {
    return(logfact(b) - logfact(a) - logfact(b - a))
  }
}

#' @title logfact function from Solow & Beet 2014
#'
#' @description
#' Helper function. Adapted from code (woodpecker-functions-version2.r)
#' provided by Adam Butler.
#'
#' @param x a non-negative integer.
#'
#' @returns a `numeric`.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @noRd

logfact <- function(x) {
  return((x > 0) * sum(log(seq(1:x))))
}
