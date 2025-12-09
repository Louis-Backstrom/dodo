#' @title Caley & Barry's (2014) "Non-constant" model
#'
#' @description
#' Non-constant population model from Caley & Barry 2014. Estimates a posterior
#' probability that the species is extant at the test time, and a point
#' estimate and one-sided \eqn{1 - \alpha} credible interval on the time of
#' extinction.
#'
#' @param records sighting records in `cbin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to the end of the observation period.
#' @param burn.in number of initial iterations to discard as burn-in (defaults
#' to 10,000).
#' @param n.iter number of iterations to run (defaults to 110,000).
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
#' @seealso [CB14B1()]
#'
#' @examples
#' # Run the fox analysis from Caley & Barry 2014
#' CB14B2(fox, init.time = 2001, test.time = 2012, n.iter = 11e4, burn.in = 1e4)
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' CB14B2(curlew$cbin, init.time = 1817, n.iter = 11e3, burn.in = 1e3)
#' }
#'
#' @export

CB14B2 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   burn.in = 10000, n.iter = 110000) {
  # Run sampler then discard burn in period
  res <- fit.func2(
    iter = n.iter, y = records, pgr.init = 0.0,
    delta.init = 0.69, eps0.init = 0, eps1.init = 1, N0.init = 1
  )
  res <- res[-(1:burn.in), ]

  # Calculate p.extant
  p.extant <- mean(is.na(res[, 6]))

  # Calculate estimate
  estimate <- median(replace(res[, 6], is.na(res[, 6]), Inf)) + init.time - 1

  # Calculate credible interval bounds
  cred.int <- as.numeric(quantile(
    replace(res[, 6], is.na(res[, 6]), Inf),
    c(0, 1 - alpha)
  )) + init.time - 1

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.extant = p.extant,
    burn.in = burn.in,
    n.iter = n.iter,
    estimate = estimate,
    cred.int = cred.int
  )

  return(output)
}
