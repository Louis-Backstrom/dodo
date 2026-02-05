#' @title Kodikara et al.'s (2020) "Uncertain" model
#'
#' @description
#' Model 2 from Kodikara et al. 2020. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-sided
#' \eqn{1 - \alpha} credible interval on the time of extinction.
#'
#' @param records sighting records in `ubin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the end of the observation period.
#' @param n.chains number of MCMC chains to run. Defaults to 4.
#' @param n.iter number of iterations in each chain. Defaults to 110,000.
#' @param n.burnin number of iterations to discard as burn-in. Defaults to
#' 10,000.
#' @param n.thin thinning rate. Defaults to 10.
#'
#' @returns a `list` object with the original parameters and the p(extant),
#' point estimate, and credible interval included as elements. The credible
#' interval is a two-element numeric vector called `cred.int`.
#'
#' @note
#' Sampling effort is assumed to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Kodikara, S., Demirhan, H., Wang, Y., Solow, A., & Stone, L. (2020).
#' Inferring extinction year using a Bayesian approach. *Methods in Ecology*
#' *and Evolution*, 11(8), 964-973. \doi{10.1111/2041-210x.13408}
#'
#' @seealso [KO20B1()]
#'
#' @examples
#' \dontrun{
#' # Run the Ivory-billed Woodpecker analysis from Kodikara et al. 2020
#' KO20B2(woodpecker$ubin, init.time = 1897, test.time = 2010)
#' # Run an example analysis using the Slender-billed Curlew data
#' KO20B2(curlew$ubin, init.time = 1817, test.time = 2022)
#' }
#'
#' @export

KO20B2 <- function(records, alpha = 0.05, init.time,
                   test.time = init.time + nrow(records) - 1, n.chains = 4,
                   n.iter = 11e4, n.burnin = 1e4, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required but could not be found!")
  }

  # Run MCMC function from Kodikara et al. 2020
  invisible(capture.output({
    posterior <- coda::mcmc.list(coda::mcmc.list(
      posterior_cer_uncer_mcmc(records[, 1], records[, 2],
        n.chains = n.chains,
        n.iter = n.iter, n.burnin = n.burnin,
        n.thin = n.thin
      )
    ))
  }))

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(posterior))
  posterior$year <- posterior$tau + init.time - 1

  # Calculate p(extant)
  p.extant <- mean(posterior$year > test.time)

  # Calculate point estimate
  estimate <- median(posterior$year)

  # Calculate credible interval bounds
  cred.int.lower <- as.numeric(quantile(posterior$year, 0))
  cred.int.upper <- as.numeric(quantile(posterior$year, 1 - alpha))

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(output)
}

# Declare certain as a known global variable (column in records data.frame):
utils::globalVariables("certain")

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

posterior_cer_uncer_mcmc <- function(y_c, y_u, n.chains = 4, n.iter = 11e4,
                                     n.burnin = 1e4, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required but could not be found!")
  }

  set.seed(1234)

  n_c <- sum(y_c)
  n_u <- sum(y_u)
  Tt <- length(y_c)

  t_n <- 0
  i <- 1
  # while (sum(y_c[i:Tt]) > 0) {
  #   t_n <- i
  #   i <- i + 1
  # }
  t_n <- max(which(y_c == 1))

  n_u_tau <- c()

  for (i in 1:Tt) {
    n_u_tau[i] <- sum(y_u[1:i])
  }

  lik <- c()

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

  model_file <- tempfile(fileext = ".txt")
  writeLines(modelStringm2, con = model_file)

  jagsModelm2 <- rjags::jags.model(
    file = model_file, data = dataList,
    n.chains = n.chains, n.adapt = n.burnin
  )
  update(jagsModelm2, n.iter = n.burnin)
  codaSamplesm2 <- rjags::coda.samples(jagsModelm2, variable.names = c(
    "tau", "pc", "pui", "puv", "theta"
  ), n.iter = n.iter, thin = n.thin)

  unlink(model_file)

  return(coda::mcmc(codaSamplesm2))
}
