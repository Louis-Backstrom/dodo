#' @title Constant Population Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S1.R).
#'
#' @param y observation data vector; 1 = observed that time period, 0 = not
#' observed
#' @param phi1 beta parameter for extinction probablity prior
#' @param phi2 beta parameter for extinction probablity prior
#' @param lambda1 beta parameter for detection probability prior
#' @param lambda2 beta parameter for detection probability prior
#' @param iter number of sampling iterations (defaults to 100,000)
#' @param initphi initial value for yearly survival probability (defaults to
#' 0.5)
#' @param initlambda initial value for yearly detection probability (given
#' extant; defaults to 0.5)
#'
#' @returns a `matrix` object containing Gibbs samples for phi, lambda and
#' z (extinction year)
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

fit.func1 <- function(y, phi1, phi2, lambda1, lambda2, iter = 1E5,
                      initphi = 0.5, initlambda = 0.5) {
  pos <- (1:length(y))[y == 1]
  final <- pos[length(pos)]
  obs <- length(y)

  curr.phi <- initphi
  curr.lambda <- initlambda

  periods <- obs - final
  result <- matrix(0, nrow = iter, ncol = 3)

  for(i in 1:iter) {
    newz <- impute(curr.phi, curr.lambda, periods)
    index <- (1:length(newz))[newz == 1]
    extra.y <- NA

    if (index > periods)  {
      extra.z <- rgeom(1, curr.phi) + 1
      curr.z <- obs + extra.z
      extra.y <- rbinom(extra.z - 1, 1, curr.lambda)
    } else {
      curr.z <- final + index
    }

    curr.phi <- rbeta(1, phi1 + 1, phi2 + curr.z - 1)
    curr.lambda <- rbeta(1, lambda1 + sum(c(y, extra.y), na.rm = TRUE),
                         lambda2 + curr.z - sum(c(y, extra.y), na.rm = TRUE))
    result[i,] <- c(curr.phi, curr.lambda, curr.z)
  }

  return(result)
}

#' @title Impute Population Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S1.R).
#'
#' @param phi extinction probability
#' @param lambda detection probability (given extant)
#' @param periods observations since last detection
#'
#' @returns a vector giving "position" of extinction since last sighting.
#' The last position is right-censored.
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

impute <- function(phi, lambda, periods) {
  result <- 1:(periods + 1)

  for(i in 1:periods) {
    result[i] <- phi * ((1 - phi) * (1 - lambda)) ^ (i - 1)
  }

  result[periods + 1] <- exp(periods * log((1 - phi) * (1 - lambda)))

  result <- result / sum(result)

  answer <- rmultinom(1, 1, result)

  return(answer)
}

#' @title sim.N Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param pgr exponential population growth rate per year
#' @param time time to end of observation period
#' @param N0 initial population size
#'
#' @returns the population trajectory from start of observation period
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

sim.N <- function(pgr, time, N0) {
  ans <- numeric(time)
  ans[1] <- N0

  if (time > 1) {
    for (t in 2:time){
      ans[t] <- ans[t - 1] * exp(pgr)
    }
  }

  return(ans)

}

#' @title Logit Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns the logit of the provided number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

logit <- function(x) {

  return(log(x / (1 - x)))

}

#' @title Inverse Logit Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns the inverse logit of the provided number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

ilogit <- function(x) {

  return(1 / (1 + exp(-x)))

}

#' @title lam.gen Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param delta per capita detection rate
#' @param N population size
#'
#' @returns a vector of yearly detection probability based on population size
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lam.gen <- function(delta, N) {

  return(1 - exp(-delta * N))

}

#' @title phi.gen Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param eps0 intercept for logit of extinction probability
#' @param eps1 population coefficient for logit of extinction probability
#' @param population size
#'
#' @returns a vector of extinction probability based on population size
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

phi.gen <- function(eps0, eps1, N) {
  lp <- eps0 - eps1 * log(N)
  if (lp > 50) {
    ans <- 1
  } else {
    ans <- ilogit(lp)
  }

  return(ans)

}

#' @title correct.negative Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a numeric vector
#'
#' @returns a numeric vector
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

correct.negative <- function(x) {
  x[sign(x) < 0] <- 0

  return(x)

}

#' @title impute.TE Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param p.cease a vector of yearly estimates of extinction during the
#' observation window
#'
#' @returns the estimated time to extinction (TE) within the observation
#' window, or otherwise returns NA if right-censored
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

impute.TE <- function(p.cease) {
  p.cease <- correct.negative(p.cease)
  p.cens <- tail(p.cease, 1)
  if (sign(p.cens < 0)) {
    p.cens <- 0
  }

  TE.draw <- rmultinom(1, 1, p.cease)

  result <- (1:length(p.cease))[TE.draw == 1]
  if (result <= (length(p.cease) - 1)) {
    return(result)
  } else {
    return(NA)
  }

}

#' @title calc.pars Function (Non-constant Population) from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param p a vector of parameters `c(pgr, delta, eps0, eps1, N0)`
#' @param y observation data vector; 1 = observed that time period, 0 = not
#' observed
#'
#' @returns a `list` object with various output parameters included
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

calc.pars <- function(p, y) {
  pgr <-  p[1]
  delta <- p[2]
  eps0 <- p[3]
  eps1 <- p[4]
  N0 <- p[5]
  obs <- length(y)

  pos <- (1:length(y))[y == 1]
  final <- tail(pos, 1)
  zero.count <- obs - final
  N.traj <- sim.N(N0 = N0, pgr = pgr, time = obs)
  lambdas <- sapply(N.traj, function(x) lam.gen(delta = delta, N = x))
  phis <- sapply(N.traj, function(x) phi.gen(eps0 = eps0, eps1 = eps1, N = x))

  p.cease <- numeric(obs + 1)
  p.cease[1:final] <- 0

  p.cease[final + 1] <- phis[final + 1]
  if (zero.count > 1) {
    for (d in (final + 2):obs) {
      log.ans <- sum(log(1 - phis[(final + 1):(d - 1)])) +
        sum(log(1 - lambdas[(final + 1):(d - 1)])) + log(phis[d])
      p.cease[d] <- exp(log.ans)
    }
  }
  p.cease <- correct.negative(p.cease)

  if (sum(p.cease) > 0) {
    log.ans <- sum(log(1 - phis[(final + 1):obs])) +
      sum(log(1 - lambdas[(final + 1):obs]))
    p.cease[obs + 1] <- exp(log.ans)
  } else {
    p.cease[obs + 1] <- 1
  }
  p.cease <- p.cease[(final + 1):(final + zero.count + 1)]
  p.cease <- p.cease / sum(p.cease)
  p.cens <- tail(p.cease, 1)

  TE <- final + impute.TE(p.cease)

  if (!is.na(TE)) {
    phis <- head(phis, TE)
    lambdas <- head(lambdas, TE)
    y.obs <- head(y, TE)
    N.traj <- head(N.traj, TE)
  } else {
    y.obs <- y
  }

  return(
    list(
      N = N.traj,
      p.cease = p.cease,
      p.cens = p.cens,
      TE = TE,
      y.obs = y.obs,
      phis = phis,
      lambdas = lambdas))

}

#' @title calc.pars.given.TE Function (Non-constant Population) from Caley &
#' Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param p a vector of parameters `c(pgr, delta, eps0, eps1, N0)`
#' @param TE a number
#' @param y observation data vector; 1 = observed that time period, 0 = not
#' observed
#'
#' @returns a `list` object with various output parameters included
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

calc.pars.given.TE <- function(p, TE, y) {
  pgr <-  p[1]
  delta <- p[2]
  eps0 <- p[3]
  eps1 <- p[4]
  N0 <- p[5]

  if(!is.na(TE)) {
    N.traj <- sim.N(N0 = N0, pgr = pgr, time = TE)
  } else {
    N.traj <- sim.N(N0 = N0, pgr = pgr, time = length(y))
  }

  lambdas <- sapply(N.traj, function(x) lam.gen(delta = delta, N = x))
  phis <- sapply(N.traj, function(x) phi.gen(eps0 = eps0, eps1 = eps1, N = x))

  return(
    list(
      N = N.traj,
      TE = TE,
      phis = phis,
      lambdas = lambdas))

}

#' @title Log-likelihood Function for lambda from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param lams a numeric vector
#' @param ys a numeric vector
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lnL.lam <- function(lams, ys){

  return(sum(log(lams[ys != 0])) + sum(log(1 - lams[ys != 1])))

}

#' @title Log-likelihood Function for phi from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param phis a numeric vector
#' @param TE a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lnL.phi <- function(phis, TE) {
  if(!is.na(TE)) {
    ans <- sum(log(1 - phis[-length(phis)])) + log(tail(phis, 1))
  } else {
    ans <- sum(log(1 - phis))
  }

  return(ans)

}

N0.prior <- function(x) {

  return(dunif(x, 5, 50, log = FALSE))

}

#' @title Log N0 Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lnL.N0.prior <- function(x) {

  return(dunif(x, 5, 50, log = TRUE))

}

#' @title pgr Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

pgr.prior <- function(x) {

  return(dunif(x, min = -2.3, max = 0.69, log = F))

}

#' @title Log pgr Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lnL.pgr.prior <- function(x) {

  return(dunif(x, min = -2.3, max = 0.69, log = TRUE))

}

#' @title delta Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

delta.prior <- function(x) {

  return(dunif(x, 0.01, 4.6, log = FALSE))

}

#' @title Log delta Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lnL.delta.prior <- function(x) {

  return(dunif(x, 0.01, 4.6, log = TRUE))

}

#' @title eps0 Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

eps0.prior <- function(x) {

  return(dunif(x, -20, 20, log = FALSE))

}

#' @title Log eps0 Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lnL.eps0.prior <- function(x) {

  return(dunif(x, -20, 20, log = TRUE))

}

#' @title eps1 Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

eps1.prior <- function(x) {

  return(dunif(x, 0, 20, log = FALSE))

}

#' @title Log eps1 Prior Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

lnL.eps1.prior <- function(x) {

  return(dunif(x, 0, 20, log = TRUE))

}

#' @title propose.N0 Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

#' @title propose.N0 Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

propose.N0 <- function(x){

  return(x)

}

#' @title propose.pgr Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#' @param sd a number (defaults to 0.05)
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

propose.pgr <- function(x, sd = 0.05) {

  return(rnorm(1, x, sd = sd))

}

#' @title propose.delta Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#' @param sd a number (defaults to 0.05)
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

propose.delta <- function(x, sd = 0.5) {

  return(rlnorm(1, log(x), sd))

}

#' @title q.delta Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x1 a number
#' @param x2 a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

q.delta <- function(x1, x2) {

  return(dlnorm(x1, log(x2), sd = 0.5, log = TRUE))

}

#' @title propose.eps0 Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#' @param sd a number (defaults to 1.5)
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

propose.eps0 <- function(x, sd = 1.5) {

  return(rnorm(1, x, sd = sd))

}

#' @title propose.eps1 Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x a number
#' @param sd a number (defaults to 0.15)
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

propose.eps1 <- function(x, sd = 0.15) {

  return(rlnorm(1, log(x), sd = sd))

}

#' @title q.eps1 Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param x1 a number
#' @param x2 a number
#'
#' @returns a number
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

q.eps1 <- function(x1, x2) {

  return(dlnorm(x1, log(x2), sd = 0.15, log = TRUE))

}

#' @title Non-constant Population Function from Caley & Barry 2014
#'
#' @description
#' Helper function. From provided code (Code S2.R).
#'
#' @param N0.init starting value for N0
#' @param y observation data vector; 1 = observed that time period, 0 = not
#' observed
#' @param iter number of iterations
#' @param pgr.init starting value for pgr
#' @param delta.init starting value for delta
#' @param eps0.init starting value for eps0
#' @param eps1.init starting value for eps1
#'
#' @returns a `matrix` object containing samples for pgr, delta, eps0, eps1,
#' N0, TE, and acceptance indicators (x4)
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @noRd

fit.func2 <- function(N0.init = 1, y, iter = 100, pgr.init = 0.0,
                      delta.init = NA, eps0.init = 0.0, eps1.init = -0.1) {
  delta.init <- -log(1 - 0.5) / N0.init

  curr.p <- c(pgr.init, delta.init, eps0.init, eps1.init, N0.init)

  result <- matrix(0 , nrow = iter, ncol = 10)
  for (i in 1:iter) {
    # cat("Doing interation", i, "of", iter, "\n")
    curr.vals <- calc.pars(p = curr.p, y = y)
    TE.imp <- curr.vals$TE
    y.obs <- curr.vals$y.obs

    curr.lnLik.delta <- lnL.lam(lams = curr.vals$lambdas, ys = y.obs) +
      lnL.delta.prior(curr.p[2])

    prop.p <- curr.p
    prop.p[2] <- propose.delta(curr.p[2])
    prop.vals <- calc.pars.given.TE(p = prop.p, TE = TE.imp, y = y)

    prop.lnLik.delta <- lnL.lam(lams = prop.vals$lambdas, ys = y.obs) +
      lnL.delta.prior(prop.p[2])

    LR <- exp(prop.lnLik.delta - curr.lnLik.delta +
                q.delta(curr.p[2], prop.p[2]) - q.delta(prop.p[2], curr.p[2]))

    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i, 8] <- 1
    } else {
      result[i, 8] <- 0
    }

    curr.lnLik.eps <- lnL.phi(phis = curr.vals$phis, TE = TE.imp) +
      lnL.eps0.prior(curr.p[3]) + lnL.eps1.prior(curr.p[4])

    prop.p <- curr.p
    prop.p[3:4] <- c(propose.eps0(curr.p[3]), propose.eps1(curr.p[4]))
    prop.vals <- calc.pars.given.TE(p = prop.p, TE = TE.imp, y = y)

    prop.lnLik.eps <- lnL.phi(phis = prop.vals$phis, TE = TE.imp) +
      lnL.eps0.prior(prop.p[3]) + lnL.eps1.prior(prop.p[4])
    LR <- exp(prop.lnLik.eps - curr.lnLik.eps + q.eps1(curr.p[4], prop.p[4]) -
                q.eps1(prop.p[4], curr.p[4]))

    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i, 9] <- 1
    } else {
      result[i, 9] <- 0
    }

    curr.lnLik.pgr <- lnL.lam(lams = curr.vals$lambdas, ys = y.obs) +
      lnL.phi(phis = curr.vals$phis, TE = TE.imp) + lnL.pgr.prior(curr.p[1])
    prop.p <- curr.p
    prop.p[1] <- propose.pgr(curr.p[1])
    prop.vals <- calc.pars.given.TE(p = prop.p, TE = TE.imp, y = y)

    prop.lnLik.pgr <- lnL.lam(lams = prop.vals$lambdas, ys = y.obs) +
      lnL.phi(phis = prop.vals$phis, TE = TE.imp) + lnL.pgr.prior(prop.p[1])
    LR <- exp(prop.lnLik.pgr - curr.lnLik.pgr)

    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i, 7] <- 1
    } else {
      result[i, 7] <- 0
    }

    curr.lnLik.N0 <- lnL.lam(lams = curr.vals$lambdas, ys = y.obs) +
      lnL.phi(phis = curr.vals$phis, TE = TE.imp) + lnL.N0.prior(curr.p[5])

    prop.p <- curr.p
    prop.p[5] <- propose.N0(curr.p[5])
    prop.vals <- calc.pars.given.TE(p = prop.p, TE = TE.imp, y = y)

    prop.lnLik.N0 <- lnL.lam(lams = prop.vals$lambdas, ys = y.obs) +
      lnL.phi(phis = prop.vals$phis, TE = TE.imp) + lnL.N0.prior(prop.p[5])

    LR <- exp(prop.lnLik.N0 - curr.lnLik.N0)

    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i, 10] <- 1
    } else {
      result[i, 10] <- 0
    }

    result[i, 1:6] <- c(curr.p,TE.imp)
  }

  return(result)

}
