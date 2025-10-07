#' @title Caley & Barry's (2014) "Constant" model
#'
#' @description
#' Constant population model from Caley & Barry 2014. Estimates a posterior
#' probability that the species is extant at the test time, and a point
#' estimate and one-tailed \eqn{1 - \alpha} credibile interval on the time of
#' extinction.
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
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
#'
#' @seealso [CB14B2()]
#'
#' @examples
#' # Run the fox analysis from Caley & Barry 2014
#' CB14B1(fox, init.time = 2001, test.time = 2012)
#'
#' @export

CB14B1 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {
  # Sort records
  records <- sort(records)

  # Create 0/1 sighting vector
  sightings <- vector(length = test.time - init.time + 1)
  sightings[records - init.time + 1] <- 1
  # if (init.time == min(records)) { # Not used, per paper
  #   sightings <- sightings[-1]
  # }

  # Set hyper priors for lambda and phi (wp for "weak prior")
  lam.wp <- c(1.1, 1.1)
  phi.wp <- c(1.1, 1.1)

  # Run sampler
  res.pp.wp <- fit.func1(
    y = sightings, phi1 = phi.wp[1], phi2 = phi.wp[2],
    lambda1 = lam.wp[1], lambda2 = lam.wp[2], iter = 1e5
  )

  # Extract extinction time posterior
  posterior <- res.pp.wp[, 3] + init.time - 1

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.extant = mean(posterior > test.time),
    estimate = median(posterior),
    cred.int = as.numeric(quantile(posterior, c(0, 1 - alpha)))
  )

  return(output)
}
