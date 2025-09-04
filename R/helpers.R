#' @title F(x) from Solow 2005
#'
#' @description
#' Helper function. Modified Equation 7 from Solow 2005.
#'
#' @param x generally \eqn{t_n} or \eqn{T}.
#' @param s sum of relative sighting times.
#' @param n number of sighting records.
#'
#' @returns a number.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#' *Mathematical Biosciences*, 195(1), 47-55. \doi{10.1016/j.mbs.2005.02.001}
#'
#' @noRd

Fx <- function(x, s, n) {
  is <- 1:floor(s / x)

  part1 <- (-1) ^ (is - 1)
  part2 <- choose(n, is)
  part3 <- (1 - (is * x / s)) ^ (n - 1)

  result <- 1 - sum(part1 * part2 * part3)
  rm(is, part1, part2, part3)

  return(result)

}

#' @title Falling factorial
#'
#' @description
#' Calculates the falling factorial function.
#'
#' @param x an integer.
#' @param j an integer.
#'
#' @returns a number.
#'
#' @noRd

ffactorial <- function(x, j) {
  return(factorial(x) / factorial(x - j))
}

#' @title Model 1 function from Kodikara et al. 2020
#'
#' @description
#' Helper function. From provided code.
#'
#' @param y_c vector of 0/1 sighting dates.
#'
#' @returns an `mcmc` object.
#'
#' @references
#' **Key Reference**
#'
#' Kodikara, S., Demirhan, H., Wang, Y., Solow, A., & Stone, L. (2020).
#' Inferring extinction year using a Bayesian approach. *Methods in Ecology*
#' *and Evolution*, 11(8), 964-973. \doi{10.1111/2041-210x.13408}
#'
#' @noRd

posterior_cer_mcmc <- function(y_c){
  set.seed(1234)

  Tt <- length(y_c)
  n <- sum(y_c)
  t_n <- 0
  i <- 1
  while(sum(y_c[i:Tt]) > 0){
    t_n <- i
    i <- i+1
  }

  lik <- c()

  dataList <- list(
    t_n = t_n,
    Tt = Tt,
    n = n,
    lik = lik
  )

  modelStringm1 <- paste0(
    "
    data {

      C <- 1000000000
      ones <- 1

    }

    model {

      for (t in 1:t_n) {
        lik[t] <- 1e-100
      }

      for (t in (t_n + 1):Tt) {
        lik[t] <- (p ^ n) * (1 - p) ^ (t - 1 - n)
      }

      lik[Tt + 1] <- (p ^ n) * (1 - p) ^ (Tt - n)

      x <- step(Tt - tau) * tau + step(tau - Tt - 1) * (Tt + 1)
      likelihood <- lik[x]

      spy <- likelihood / C
      ones ~ dbern(spy)

      tau_0 ~ dnegbin(theta, 1)
      tau <- tau_0 + 1
      theta ~ dunif(0, 1)
      p ~ dunif(0, 1)

    }
    "
    )

  writeLines(modelStringm1, con = "model_m1.txt")

  jagsModelm1 <- rjags::jags.model(file = "model_m1.txt", data = dataList,
                                   n.chains = 4, n.adapt = 60000)
  update(jagsModelm1, n.iter = 60000)
  codaSamplesm1 <- rjags::coda.samples(jagsModelm1, variable.names = c(
    "tau", "p", "theta"), n.iter = 130000, thin = 13)
  unlink("model_m1.txt")

  return(coda::mcmc(codaSamplesm1))

}

#' @title Model 2 function from Kodikara et al. 2020
#'
#' @description
#' Helper function. From provided code.
#'
#' @param y_c vector of 0/1 certain sighting dates.
#' @param y_u vector of 0/1 uncertain sighting dates
#'
#' @returns an `mcmc` object.
#'
#' @references
#' **Key Reference**
#'
#' Kodikara, S., Demirhan, H., Wang, Y., Solow, A., & Stone, L. (2020).
#' Inferring extinction year using a Bayesian approach. *Methods in Ecology*
#' *and Evolution*, 11(8), 964-973. \doi{10.1111/2041-210x.13408}
#'
#' @noRd

posterior_cer_uncer_mcmc <- function(y_c,y_u){
  set.seed(1234)

  n_c <- sum(y_c)
  n_u <- sum(y_u)
  Tt <- length(y_c)

  t_n <- 0
  i <- 1
  while(sum(y_c[i:Tt]) > 0){
    t_n <- i
    i <- i+1
  }

  n_u_tau <- c()

  for(i in 1:Tt){
    n_u_tau[i] <- sum(y_u[1:i])
  }

  lik<-c()

  dataList <- list(
    t_n = t_n,
    Tt = Tt,
    nc = n_c,
    n_u = n_u,
    n_u_tau = n_u_tau,
    lik = lik
  )

  modelStringm2 <- paste0(
    "
    data {

      C <- 1000000000
      ones <- 1

    }

    model {

      for (t in 1:t_n) {
        lik[t] <- 1e-100
      }

      for (t in (t_n + 1):Tt) {
        lik[t] <-
          (pc ^ nc) *
          (((1 - pc) * pu) ^ n_u_tau[t - 1]) *
          (((1 - pc) * (1 - pu)) ^ (t - 1 - nc - n_u_tau[t - 1])) *
          (pui ^ (n_u - n_u_tau[t - 1])) *
          ((1 - pui) ^ (Tt - (t - 1) - (n_u - n_u_tau[t - 1])))
      }

      lik[Tt + 1] <-
        (pc ^ nc) *
        (((1 - pc) * pu) ^ n_u) *
        (((1 - pc) * (1 - pu)) ^ (Tt - nc - n_u))

      x <- step(Tt - tau) * tau + step(tau - Tt - 1) * (Tt + 1)
      likelihood <- lik[x]

      spy <- likelihood / C
      ones ~ dbern(spy)

      tau_0 ~ dnegbin(theta, 1)
      tau <- tau_0 + 1
      theta ~ dunif(0, 1)

      pu <- puv * (1 - pui) + pui * (1 - puv) + puv * pui
      pui ~ dunif(0, 1)
      puv ~ dunif(0, 1)
      pc ~ dunif(0, 1)

    }
    "
    )

  writeLines(modelStringm2, con = "model_m2.txt")

  jagsModelm2 <- rjags::jags.model(file = "model_m2.txt", data = dataList,
                                   n.chains = 4, n.adapt = 60000)
  update(jagsModelm2, n.iter = 60000)
  codaSamplesm2 <- rjags::coda.samples(jagsModelm2, variable.names = c(
    "tau", "pc", "pui", "puv", "theta"), n.iter = 130000, thin = 13)
  unlink("model_m2.txt")
  return(coda::mcmc(codaSamplesm2))

}

