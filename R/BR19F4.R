#' @title Brook et al.'s (2019) "McInerny" model
#'
#' @description
#' The McInerny et al. 2006 variant of the model from Brook et al. 2019.
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
#' McInerny, G. J., Roberts, D. L., Davy, A. J., & Cribb, P. J. (2006).
#' Significance of sighting rate in inferring extinction and threat.
#' *Conservation Biology*, 20(2), 562-567.
#' \doi{10.1111/j.1523-1739.2006.00377.x}
#'
#' @seealso [BR19F1()], [BR19F2()], [BR19F3()], [BR19F5()]
#'
#' @examples
#' # Run the "Extreme" Ivory-billed Woodpecker analysis from Brook et al. 2019
#' BR19F4(woodpecker$ucon, alpha = 0.1, test.time = 2010, n.iter = 1e3)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' BR19F4(curlew$ucon, test.time = 2022, n.iter = 1e3)
#' }
#'
#' @export

BR19F4 <- function(records, alpha = 0.05, init.time = min(records$time),
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
    ]
  )

  # Run MC06F1 on samples
  estimates <- sapply(samples, function(x) {
    MC06F1(
      records = init.time:test.time %in% x, alpha = 0.5,
      init.time = init.time
    )$conf.int[2]
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
