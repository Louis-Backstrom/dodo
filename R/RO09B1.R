#' @title Rout et al.'s (2009) "Declining" model
#'
#' @description
#' Equation 6 and others from Rout et al. 2009. Estimates a pre-extinction rate
#' of decline, a Bayes factor comparing competing hypotheses of extinction
#' / persistence, and a posterior probability that the species is extant at the
#' end of the observation period.
#'
#' @param records sighting records in `cbin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param pi prior probability that \eqn{H_0} is true (defaults to
#' \eqn{\pi = 0.5}).
#' @param n.chains number of MCMC chains to run. Defaults to 4.
#' @param n.iter number of iterations in each chain. Defaults to 110,000.
#' @param n.burnin number of iterations to discard as burn-in. Defaults to
#' 10,000.
#'
#' @returns a `list` object with the original parameters and the rate of
#' decline (\eqn{a}), the Bayes factor, and p(extant) included as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant. The Bayes Factor presented here is the inverse of the Bayes
#' Factor as presented in the original paper, to allow for comparability with
#' other models in this package (values > 1 imply extinction).
#'
#' @references
#' **Key Reference**
#'
#' Rout, T. M., Salomon, Y., & McCarthy, M. A. (2009). Using sighting records
#' to declare eradication of an invasive species. *Journal of Applied Ecology*,
#' 46(1), 110-117. \doi{10.1111/j.1365-2664.2008.01586.x}
#'
#' **Other References**
#'
#' Solow, A. R. (1993). Inferring Extinction from Sighting Data. *Ecology*,
#' 74(3), 962-964. \doi{10.2307/1940821}
#'
#' @examples
#' \dontrun{
#' # Run the bitterweed analysis from Rout et al. 2009
#' RO09B1(bitterweed)
#' # Run an example analysis using the Slender-billed Curlew data
#' RO09B1(curlew$cbin)
#' }
#'
#' @export

RO09B1 <- function(records, pi = 0.5, n.chains = 4, n.iter = 11e4,
                   n.burnin = 1e4) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required but could not be found!")
  }

  # Set up JAGS model
  data <- list(
    records = records,
    N = length(records)
  )

  model_string <- "
    model {
      a ~ dunif(0, 1)
      m ~ dnorm(0, 1.0E-6)
      prec ~ dgamma(0.001, 0.001)

      for (i in 1:N) {
        meanlam[i] <- m * pow(i, -a)
        records[i] ~ dnorm(meanlam[i], prec)
      }
    }
    "

  inits <- rep(list(list(a = 0.5, m = 0.5, prec = 100)), n.chains)

  # Run the JAGS model
  model <- rjags::jags.model(textConnection(model_string),
    data = data,
    inits = inits,
    n.chains = n.chains,
    n.adapt = n.burnin,
    quiet = TRUE
  )

  # Sink (to suppress hyper-verbose console outputs)
  sink(file = tempfile())

  update(model, n.burnin)

  samples <- rjags::coda.samples(model,
    variable.names = c("a", "m", "prec"),
    n.iter = n.iter
  )

  on.exit(sink(), add = TRUE)

  a_mean <- summary(samples)$statistics["a", "Mean"]
  S <- length(records)
  n <- sum(records)
  sn <- max(which(records == 1))

  # Calculate Bayes factor
  Bayes.factor <- (n * (1 - a_mean) - 1) /
    ((S / sn)^(n * (1 - a_mean) - 1) - 1)

  # Calculate posterior
  posterior <- (1 + ((1 - pi) / (pi * Bayes.factor)))^-1

  # Output
  output <- list(
    records = records,
    pi = pi,
    a = summary(samples)$statistics["a", "Mean"],
    Bayes.factor = 1 / Bayes.factor,
    p.extant = posterior
  )

  return(output)
}
