#' @title Wang et al.'s (2016) "Adaptive" model
#'
#' @description
#' The Adaptive Beta Method from Wang et al. 2016. Estimates a posterior
#' distribution on the time of extinction, with associated point estimate and
#' one-sided \eqn{1 - \alpha} credible interval.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param use.mpfr whether or not to use multiple-precision floating-point
#' computation via the `Rmpfr` package. Defaults to `FALSE`, but is necessary
#' for larger datasets (N > 50). Note that this is very slow, and experimental!
#'
#' @returns a `list` object with the original parameters and the point estimate
#' and credible interval included as elements. The credible interval is a
#' two-element numeric vector called `cred.int`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' @examples
#' # Run the Anabarella analysis from Wang et al. 2016
#' WA16B1(-anabarella, alpha = 0.1, use.mpfr = FALSE)
#' \dontrun{
#' WA16B1(-anabarella, alpha = 0.1, use.mpfr = TRUE) # Slow!
#' }
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' WA16B1(curlew$ccon, init.time = 1817, use.mpfr = FALSE) # Fails!
#' }
#' \dontrun{
#' WA16B1(curlew$ccon, init.time = 1817, use.mpfr = TRUE) # Slow!
#' }
#'
#' @export

WA16B1 <- function(records, alpha = 0.05, init.time = min(records),
                   use.mpfr = FALSE) {
  # Sort records
  records <- sort(records)

  # Determine base for ABM function
  if (init.time == min(records)) {
    base <- NULL
  } else {
    base <- init.time
  }

  # Run the ABM model
  abm_results <- abm(
    x = records, distance = FALSE, ext = FALSE, base = base,
    prmean = 0, prSD = 2, alpha = alpha, PLOT = 0, use.mpfr = use.mpfr
  )

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    estimate = abm_results[1],
    cred.int = abm_results[2:3]
  )

  return(output)
}
