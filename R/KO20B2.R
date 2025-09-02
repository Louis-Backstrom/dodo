#' @title Kodikara et al.'s (2020) "Uncertain" Model
#'
#' @description
#' Model 2 from Kodikara et al. 2020. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-tailed
#' \eqn{1 - \alpha} credibile interval on the time of extinction.
#'
#' @param records `data.frame` with two columns: `time` and `certain`. The
#' `time` column contains the date of all sightings, which are each either
#' certain (`certain = TRUE`) or uncertain (`certain = FALSE`).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
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
#' @export

KO20B2 <- function(records, alpha = 0.05, init.time = min(records$time),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {

  # Sort records
  records <- sort_by(records, ~time)

  # Create 0/1 sighting vectors
  sightings_c <- vector(length = test.time - init.time + 1)
  sightings_c[subset(records, certain == TRUE)$time - init.time + 1] <- 1
  if (init.time == min(records$time)) {
    sightings_c <- sightings_c[-1]
  }
  sightings_u <- vector(length = test.time - init.time + 1)
  sightings_u[subset(records, certain == FALSE)$time - init.time + 1] <- 1
  if (init.time == min(records$time)) {
    sightings_u <- sightings_u[-1]
  }

  # Sink (to suppress hyper-verbose console outputs)
  sink(file = tempfile())

  # Run MCMC function from Kodikara et al. 2020
  posterior <- coda::mcmc.list(coda::mcmc.list(
    posterior_cer_uncer_mcmc(sightings_c, sightings_u)))

  sink()

  # Extract posteriors
  posterior <- as.data.frame(coda:::as.matrix.mcmc.list(posterior))

  # Calculate p(extant)
  p.extant <- mean(posterior$tau + init.time > test.time)

  # Calculate point estimate
  estimate <- median(posterior$tau) + init.time

  # Calculate credible interval bounds
  cred.int.lower <- as.numeric(quantile(posterior$tau, 0)) + init.time
  cred.int.upper <- as.numeric(quantile(posterior$tau, 1 - alpha)) + init.time

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.extant = p.extant,
    estimate = estimate,
    conf.int = c(cred.int.lower, cred.int.upper)
  )

  return(output)

}
