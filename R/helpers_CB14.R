fit.func1 <- function(y, phi1, phi2, lambda1, lambda2, iter = 1E5,
                      initphi = 0.5, initlambda = 0.5) {
  pos <- (1:length(y))[y == 1]
  final <- pos[length(pos)]
  obs <- length(y)

  curr.phi <- initphi
  curr.lambda <- initlambda

  periods <- obs - final
  result <- matrix(0, nrow = iter, ncol = 3)

  for (i in 1:iter) {
    newz <- impute(curr.phi, curr.lambda, periods)
    index <- (1:length(newz))[newz == 1]
    extra.y <- NA

    if (index > periods) {
      extra.z <- rgeom(1, curr.phi) + 1
      curr.z <- obs + extra.z
      extra.y <- rbinom(extra.z - 1, 1, curr.lambda)
    } else {
      curr.z <- final + index
    }

    curr.phi <- rbeta(1, phi1 + 1, phi2 + curr.z - 1)
    curr.lambda <- rbeta(
      1, lambda1 + sum(c(y, extra.y), na.rm = TRUE),
      lambda2 + curr.z - sum(c(y, extra.y), na.rm = TRUE)
    )
    result[i, ] <- c(curr.phi, curr.lambda, curr.z)
  }

  return(result)
}

impute <- function(phi, lambda, periods) {
  probs <- phi * ((1 - phi) * (1 - lambda))^(0:(periods - 1))
  probs <- c(probs, ((1 - phi) * (1 - lambda))^periods)
  probs <- probs / sum(probs)
  return(as.vector(rmultinom(1, 1, probs)))
}

sim.N <- function(pgr, time, N0) {
  return(exp(log(N0) + cumsum(c(0, rep(pgr, time - 1)))))
}

logit <- function(x) {
  return(log(x / (1 - x)))
}

correct.negative <- function(x) {
  return(pmax(x, 0))
}

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

calc.pars <- function(p, y) {
  pgr <- p[1]
  delta <- p[2]
  eps0 <- p[3]
  eps1 <- p[4]
  N0 <- p[5]
  obs <- length(y)

  pos <- (1:length(y))[y == 1]
  final <- tail(pos, 1)
  zero.count <- obs - final
  N.traj <- sim.N(N0 = N0, pgr = pgr, time = obs)
  lambdas <- 1 - exp(-delta * N.traj)
  lp <- eps0 - eps1 * log(N.traj)
  phis <- ifelse(lp > 50, 1, plogis(lp))

  p.cease <- numeric(obs + 1)
  p.cease[1:final] <- 0

  if (zero.count >= 1) {
    phis_sub <- phis[(final + 1):(final + zero.count)]
    lambdas_sub <- lambdas[(final + 1):(final + zero.count)]

    if (zero.count == 1) {
      p_cease_vec <- phis_sub
    } else {
      cp <- cumprod((1 - phis_sub[-zero.count]) * (1 - lambdas_sub[-zero.count]))
      p_cease_vec <- numeric(zero.count)
      p_cease_vec[1] <- phis_sub[1]
      p_cease_vec[2:zero.count] <- phis_sub[2:zero.count] * cp
    }

    p.cease[(final + 1):(final + zero.count)] <- p_cease_vec
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
      lambdas = lambdas
    )
  )
}

calc.pars.given.TE <- function(p, TE, y) {
  pgr <- p[1]
  delta <- p[2]
  eps0 <- p[3]
  eps1 <- p[4]
  N0 <- p[5]

  if (!is.na(TE)) {
    N.traj <- sim.N(N0 = N0, pgr = pgr, time = TE)
  } else {
    N.traj <- sim.N(N0 = N0, pgr = pgr, time = length(y))
  }
  lambdas <- 1 - exp(-delta * N.traj)
  lp <- eps0 - eps1 * log(N.traj)
  phis <- ifelse(lp > 50, 1, plogis(lp))

  return(
    list(
      N = N.traj,
      TE = TE,
      phis = phis,
      lambdas = lambdas
    )
  )
}

lnL.lam <- function(lams, ys) {
  return(sum(log(lams[ys != 0])) + sum(log(1 - lams[ys != 1])))
}

lnL.phi <- function(phis, TE) {
  if (!is.na(TE)) {
    ans <- sum(log(1 - phis[-length(phis)])) + log(tail(phis, 1))
  } else {
    ans <- sum(log(1 - phis))
  }

  return(ans)
}

N0.prior <- function(x) {
  return(dunif(x, 5, 50, log = FALSE))
}

lnL.N0.prior <- function(x) {
  return(dunif(x, 5, 50, log = TRUE))
}

pgr.prior <- function(x) {
  return(dunif(x, min = -2.3, max = 0.69, log = F))
}

lnL.pgr.prior <- function(x) {
  return(dunif(x, min = -2.3, max = 0.69, log = TRUE))
}

delta.prior <- function(x) {
  return(dunif(x, 0.01, 4.6, log = FALSE))
}

lnL.delta.prior <- function(x) {
  return(dunif(x, 0.01, 4.6, log = TRUE))
}

eps0.prior <- function(x) {
  return(dunif(x, -20, 20, log = FALSE))
}

lnL.eps0.prior <- function(x) {
  return(dunif(x, -20, 20, log = TRUE))
}

eps1.prior <- function(x) {
  return(dunif(x, 0, 20, log = FALSE))
}

lnL.eps1.prior <- function(x) {
  return(dunif(x, 0, 20, log = TRUE))
}

propose.N0 <- function(x) {
  return(x)
}

propose.pgr <- function(x, sd = 0.05) {
  return(rnorm(1, x, sd = sd))
}

propose.delta <- function(x, sd = 0.5) {
  return(rlnorm(1, log(x), sd))
}

q.delta <- function(x1, x2) {
  return(dlnorm(x1, log(x2), sdlog = 0.5, log = TRUE))
}

propose.eps0 <- function(x, sd = 1.5) {
  return(rnorm(1, x, sd = sd))
}

propose.eps1 <- function(x, sd = 0.15) {
  return(rlnorm(1, log(x), sdlog = sd))
}

q.eps1 <- function(x1, x2) {
  return(dlnorm(x1, log(x2), sdlog = 0.15, log = TRUE))
}

fit.func2 <- function(N0.init = 1, y, iter = 100, pgr.init = 0.0,
                      delta.init = NA, eps0.init = 0.0, eps1.init = -0.1) {
  delta.init <- -log(1 - 0.5) / N0.init

  curr.p <- c(pgr.init, delta.init, eps0.init, eps1.init, N0.init)

  result <- matrix(0, nrow = iter, ncol = 10)
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

    if (runif(1) < LR & !is.na(LR)) {
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

    if (runif(1) < LR & !is.na(LR)) {
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

    if (runif(1) < LR & !is.na(LR)) {
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

    if (runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i, 10] <- 1
    } else {
      result[i, 10] <- 0
    }

    result[i, 1:6] <- c(curr.p, TE.imp)
  }

  return(result)
}
