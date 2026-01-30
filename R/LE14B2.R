#' @title Lee et al.'s (2014) "Variable reliability" model
#'
#' @description
#' The model from Lee et al. 2014, incorporating uncertain sightings. Estimates
#' a posterior probability that the species is extant at the end of the
#' observation period, and a point estimate and one-tailed \eqn{1 - \alpha}
#' credible interval on the time of extinction, conditional on extinction
#' occurring before the end of the observation period.
#'
#' @param records sighting records in `ubin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param n.chains number of MCMC chains to run. Defaults to 4.
#' @param n.iter number of iterations in each chain. Defaults to 15,000.
#' @param n.burnin number of iterations to discard as burn-in. Defaults to
#' 5,000.
#' @param n.thin thinning rate. Defaults to 10.
#'
#' @returns a `list` object with the original parameters and the p(extant),
#' point estimate, and credible interval included as elements. The credible
#' interval is a two-element numeric vector called `cred.int`.
#'
#' @note
#' Sampling effort is assumed to be constant. Uses JAGS instead of BUGS.
#'
#' @references
#' **Key Reference**
#'
#' Lee, T. E., McCarthy, M. A., Wintle, B. A., Bode, M., Roberts, D. L., &
#' Burgman, M. A. (2014). Inferring extinctions from sighting records of
#' variable reliability. *Journal of Applied Ecology*, 51(1), 251-258.
#' \doi{10.1111/1365-2664.12144}
#'
#' @examples
#' \dontrun{
#' # Run the example analysis from Lee et al. 2014
#' LE14B2(lee_s1, init.time = 1)
#' # Run an example analysis using the Slender-billed Curlew data
#' LE14B2(curlew$ubin, init.time = 1817)
#' }
#'
#' @export

LE14B2 <- function(records, alpha = 0.05, init.time, n.chains = 4,
                   n.iter = 15000, n.burnin = 5000, n.thin = 10) {
  # Check if rjags is installed
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required but could not be found!")
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required but could not be found!")
  }

  # Convert records into model format
  data <- list(y = records$certain, z = records$uncertain, T = nrow(records))

  # JAGS model
  model <- "
  model {
    te ~ dunif(0, T)
    m[1] ~ dunif(0, 100)
    f[1] <- 0
    m[2] ~ dunif(0, 100)
    f[2] ~ dunif(0, 100)
    Extant ~ dbern(0.5)

    for (i in 1:T) {
      mm1[i] <- Extant * m[1] + (1 - Extant) * step(te - i) * m[1] + f[1]
      p1[i]  <- max(1.0E-6, min(1 - 1.0E-6, 1 - exp(-mm1[i])))
      y[i] ~ dbern(p1[i])

      mm2[i] <- Extant * m[2] + (1 - Extant) * step(te - i) * m[2] + f[2]
      p2[i]  <- max(1.0E-6, min(1 - 1.0E-6, 1 - exp(-mm2[i])))
      z[i] ~ dbern(p2[i])
    }
  }
  "

  # Initialize and compile model
  jags_model <- rjags::jags.model(
    textConnection(model),
    data = data,
    n.chains = n.chains,
    quiet = TRUE
  )

  # Burn-in
  update(jags_model, n.iter = n.burnin)

  # Specify parameters
  parameters <- c("te", "m", "f", "Extant")

  # Sample from posterior
  samples <- rjags::coda.samples(
    model = jags_model,
    variable.names = parameters,
    n.iter = n.iter - n.burnin,
    thin = n.thin
  )

  # Convert to matrix
  sims <- as.matrix(samples)

  # Extract posteriors (only keep te estimates for samples where Extant == 0)
  p.extant <- mean(sims[, "Extant"])
  te <- sims[sims[, "Extant"] == 0, "te"]
  estimate <- mean(te) + init.time - 1
  cred.int <- as.numeric(quantile(te, c(0, 0.95)) + init.time - 1)

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    n.chains = n.chains,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = cred.int
  )

  return(output)
}
