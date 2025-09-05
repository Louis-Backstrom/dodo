#' @title Kodikara et al.'s (2020) "Certain-only" model
#'
#' @description
#' Model 1 from Kodikara et al. 2020. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-tailed
#' \eqn{1 - \alpha} credibile interval on the time of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
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
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Kodikara, S., Demirhan, H., Wang, Y., Solow, A., & Stone, L. (2020).
#' Inferring extinction year using a Bayesian approach. *Methods in Ecology*
#' *and Evolution*, 11(8), 964-973. \doi{10.1111/2041-210x.13408}
#'
#' @seealso [KO20B2()]
#'
#' @examples
#' # Run the Ivory-billed Woodpecker analysis from Kodikara et al. 2020
#' KO20B1(woodpecker1, test.time = 2010)
#'
#' @export

KO20B1 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {

  # Sort records
  records <- sort(records)

  # Create 0/1 sighting vector
  sightings <- vector(length = test.time - init.time + 1)
  sightings[records - init.time + 1] <- 1
  if (init.time == min(records)) {
    sightings <- sightings[-1]
  }

  # Sink (to suppress hyper-verbose console outputs)
  sink(file = tempfile())

  # Run MCMC function from Kodikara et al. 2020
  posterior <- coda::mcmc.list(coda::mcmc.list(posterior_cer_mcmc(sightings)))

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
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(output)

}
