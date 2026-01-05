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
#' @param n.iter number of iterations in each chain. Defaults to 11000.
#' @param n.burnin number of iterations to discard as burn-in. Defaults to 1000.
#' @param n.thin thinning rate. Defaults to 10.
#' @param debug whether to open the OpenBUGS interface during execution.
#' Defaults to `FALSE`.
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
#' Lee, T. E., McCarthy, M. A., Wintle, B. A., Bode, M., Roberts, D. L., &
#' Burgman, M. A. (2014). Inferring extinctions from sighting records of
#' variable reliability. *Journal of Applied Ecology*, 51(1), 251-258.
#' \doi{10.1111/1365-2664.12144}
#'
#' @examples
#' # Run the example analysis from Lee et al. 2014
#' LE14B2(lee_s1, init.time = 1)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' LE14B2(curlew$ubin, init.time = 1817)
#' }
#'
#' @export

LE14B2 <- function(records, alpha = 0.05, init.time, n.chains = 4,
                   n.iter = 15000, n.burnin = 5000, n.thin = 10,
                   debug = FALSE) {
  # Check if R2OpenBUGS is installed
  if (!requireNamespace("R2OpenBUGS", quietly = TRUE)) {
    stop("Package 'R2OpenBUGS' is required but could not be found!")
  }

  # Convert records into model format
  data <- list(y = records$certain, z = records$uncertain, T = nrow(records))

  # Define model
  model <- "
  model {
    te ~ dunif(0, T)
    m[1] ~ dunif(0, 100)
    f[1] <- 0
    m[2] ~ dunif(0, 100)
    f[2] ~ dunif(0, 100)
    Extant ~ dbern(0.5)

    for (i in 1:T) {
      mm1[i] <- Extant * m[1] + (1-Extant) * step(te - i) * m[1] + f[1]
      p1[i] <- 1 - exp(-mm1[i])
      y[i] ~ dbern(p1[i])

      mm2[i] <- Extant * m[2] + (1-Extant) * step(te - i) * m[2] + f[2]
      p2[i] <- 1 - exp(-mm2[i])
      z[i] ~ dbern(p2[i])
    }
  }
  "

  model_file <- tempfile(fileext = ".txt")
  writeLines(model, model_file)

  # Specify parameters
  parameters <- c("te", "m", "f", "Extant")

  # Run the chains
  bugs_model <- R2OpenBUGS::bugs(
    data = data,
    inits = NULL,
    parameters.to.save = parameters,
    model.file = model_file,
    n.chains = n.chains,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    debug = debug
  )

  unlink(model_file)

  # Extract posteriors (only keep te estimates for samples where Extant == 0)
  p.extant <- mean(bugs_model$sims.list$Extant)
  te <- bugs_model$sims.list$te[bugs_model$sims.list$Extant == 0]

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
