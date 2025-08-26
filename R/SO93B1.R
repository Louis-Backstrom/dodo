#' @title Solow's (1993) "Bayesian" Model
#'
#' @description
#' Equation 3 and others from Solow 1993. Estimates a Bayes factor and
#' posterior distribution on probability that the species is extant, with
#' associated point estimate and two-sided credible interval.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case  this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param prior prior distribution of extinction probability.
#' Defaults to a \eqn{U(0, 1)} distribution.
#' @param keep.dists whether to include the posterior distribution vectors in
#' the final function output or not (defaults to `FALSE`).
#'
#' @returns a `list` object with the original parameters and the Bayes factor,
#' point estimate, and credible interval included as elements. The credible
#' interval is a two-element numeric vector called `cred.int`. The prior and
#' posterior distribution vectors are also included, if `keep.dists` is set
#' to `TRUE`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (1993). Inferring Extinction from Sighting Data. *Ecology*,
#' 74(3), 962-964. \doi{10.2307/1940821}
#'
#' @export

SO93B1 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   prior = runif(n = 1e6, min = 0, max = 1),
                   keep.dists = FALSE) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)
  if (init.time == min(records)) {
    n <- n - 1
  }

  # Determine length of sighting record
  tn <- max(records) - init.time

  # Determine current time
  bigT <- test.time - init.time

  # Calculate Bayes factor
  Bayes.factor <- (n - 1) / ((bigT / tn) ^ (n - 1) - 1)

  # Sample from posterior
  posterior <- (1 + ((1 - prior) / (prior * Bayes.factor))) ^ -1

  # Output

  if (keep.dists == TRUE) {
    output <- list(
      records = records,
      alpha = alpha,
      init.time = init.time,
      test.time = test.time,
      prior = prior,
      posterior = posterior,
      Bayes.factor = Bayes.factor,
      estimate = mean(posterior),
      cred.int = as.vector(quantile(posterior, c(0.025, 0.975)))
    )
  } else {
    output <- list(
      records = records,
      alpha = alpha,
      init.time = init.time,
      test.time = test.time,
      Bayes.factor = Bayes.factor,
      estimate = mean(posterior),
      cred.int = as.vector(quantile(posterior, c(0.025, 0.975)))
    )
  }

  return(output)

}
