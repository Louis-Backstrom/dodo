#' @title Brook et al.'s (2019) "Solow & Roberts" model
#'
#' @description
#' The Solow & Roberts 2003 variant of the model from Brook et al. 2019.
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
#' **Other References**
#'
#' Solow, A. R., & Roberts, D. L. (2003). A nonparametric test for extinction
#' based on a sighting record. *Ecology*, 84(5), 1329-1332.
#' \doi{10.1890/0012-9658(2003)084[1329:ANTFEB]2.0.CO;2}
#'
#' Robson, D. S., & Whitlock, J. H. (1964). Estimation of a truncation point.
#' *Biometrika*, 51(1-2), 33-39. \doi{10.1093/biomet/51.1-2.33}
#'
#' @seealso [BR19F1()], [BR19F3()], [BR19F4()], [BR19F5()]
#'
#' @examples
#' # Run the "Extreme" Ivory-billed Woodpecker analysis from Brook et al. 2019
#' BR19F2(woodpecker$ucon, alpha = 0.1, test.time = 2010, n.iter = 1e3)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' BR19F2(curlew$ucon, test.time = 2022, n.iter = 1e3)
#' }
#'
#' @export

BR19F2 <- function(records, alpha = 0.05, init.time = min(records$time),
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
      runif(nrow(records)) <= records$certainty,
      "time"
    ],
    simplify = FALSE
  )

  # Run SR03F1 on samples
  estimates <- sapply(samples, function(x) {
    SR03F1(
      records = x, alpha = alpha,
      test.time = test.time
    )$estimate
  })

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
