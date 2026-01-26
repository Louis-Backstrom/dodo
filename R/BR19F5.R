#' @title Brook et al.'s (2019) "Last Appearance Date" model
#'
#' @description
#' The Last Appearance Date variant of the model from Brook et al. 2019.
#' Estimates a "probability" of extinction at the test time, and a two-sided
#' \eqn{1 - \alpha} confidence interval and point estimate on the time of
#' extinction. Sighting uncertainty is incorporated.
#'
#' @param records sighting records in `ucon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param n.iter number of iterations to run (defaults to 10,000).
#'
#' @returns a `list` object with the original parameters and the probability,
#' point estimate, and confidence interval included as elements. The confidence
#' interval is a two-element numeric vector called `conf.int`.
#'
#' @note
#' Sampling effort is assumed to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Brook, B. W., Buettel, J. C., & Jaric, I. (2019). A fast re-sampling method
#' for using reliability ratings of sightings with extinction-date estimators.
#' *Ecology*, 100(9), e02787. \doi{10.1002/ecy.2787}
#'
#' @seealso [BR19F1()], [BR19F2()], [BR19F3()], [BR19F4()]
#'
#' @examples
#' # Run the "Extreme" Ivory-billed Woodpecker analysis from Brook et al. 2019
#' BR19F5(woodpecker$ucon, alpha = 0.1, test.time = 2010, n.iter = 1e3)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' BR19F5(curlew$ucon, test.time = 2022, n.iter = 1e3)
#' }
#'
#' @export

BR19F5 <- function(records, alpha = 0.05, init.time = min(records$time),
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   n.iter = 1e4) {
  # Sort records
  records <- sort_by(records, ~time)

  # If using first record as init.time, remove this from the record sequence
  if (init.time == min(records$time)) {
    records <- records[-1, ]
  }

  # Sample using sighting certainty as inclusion probability
  samples <- replicate(
    n = n.iter,
    expr = records[
      runif(nrow(records)) < records$certainty,
      "time"
    ]
  )

  # Run LAD estimation on samples
  estimates <- sapply(samples, max)

  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    n.iter = n.iter,
    p.extant = mean(estimates > test.time),
    estimate = median(estimates),
    conf.int = as.numeric(quantile(estimates, c(0, 1 - alpha)))
  )

  return(output)
}
