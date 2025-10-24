#' @title Solow's (1993) "Bayesian" model
#'
#' @description
#' Equation 3 and others from Solow 1993. Estimates a Bayes factor and
#' posterior distribution on probability that the species is extant, with
#' associated point estimate and two-sided \eqn{1 - \alpha} credible interval.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param pi prior probability that \eqn{H_0} is true (defaults to
#' \eqn{\pi = 0.5}).
#'
#' @returns a `list` object with the original parameters and the Bayes factor
#' and p(extant) included as elements.
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
#' @seealso [SO93F1()]
#'
#' @examples
#' # Run the Caribbean Monk Seal analysis from Solow 1993
#' SO93B1(monk_seal, test.time = 1992)
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' SO93B1(curlew$ccon, init.time = 1817, test.time = 2022)
#' }
#'
#' @export

SO93B1 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   pi = 0.5) {
  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # If using first record as init.time, remove this from the record sequence
  if (init.time == min(records)) {
    records <- tail(records, -1)
    n <- n - 1
  }

  # Determine length of sighting record
  tn <- max(records) - init.time

  # Determine current time
  bigT <- test.time - init.time

  # Calculate Bayes factor
  Bayes.factor <- (n - 1) / ((bigT / tn)^(n - 1) - 1)

  # Calculate posterior
  posterior <- (1 + ((1 - pi) / (pi * Bayes.factor)))^-1

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    pi = pi,
    Bayes.factor = Bayes.factor,
    p.extant = posterior
  )

  return(output)
}
