#' @title Brook et al.'s (2019) "Roberts & Solow" model
#'
#' @description
#' The Roberts & Solow 2003 variant of the model from Brook et al. 2019.
#' Estimates a "probability" of extinction at the test time, and a two-tailed
#' \eqn{1 - \alpha} confidence interval and point estimate on the time of
#' extinction. Sighting uncertainty is incorporated.
#'
#' @param records sighting records in `ucon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param cores number of cores to use (defaults to `NULL`, in which case
#' `parallel::detectCores()` is run).
#' @param n.iter number of iterations to run (defaults to 10000).
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
#' Roberts, D. L., & Solow, A. R. (2003). Flightless birds: when did the dodo
#' become extinct? *Nature*, 426(6964), 245. \doi{10.1038/426245a}
#'
#' @seealso [BR19F1()], [BR19F2()], [BR19F4()], [BR19F5()]
#'
#' @examples
#' # Run the "Extreme" Ivory-billed Woodpecker analysis from Brook et al. 2019
#' BR19F3(woodpecker$ucon,
#'   test.time = 2010, alpha = 0.1, cores = 2,
#'   n.iter = 1e3
#' )
#' # NB: alpha = 0.1 as this package presents two-sided confidence intervals
#' # vs. one-sided confidence intervals in the original paper; the upper bound
#' # of the one-sided alpha = 0.05 CI in Brook et al. is equivalent to the
#' # upper bound of the two-sided alpha = 0.1 CI in this example.
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' BR19F3(curlew$ucon, test.time = 2022, cores = 2, n.iter = 1e3)
#' }
#'
#' @export

BR19F3 <- function(records, alpha = 0.05,
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   cores = NULL, n.iter = 1e4) {
  # Sort records
  records <- sort_by(records, ~time)

  # Rename columns to fit model function
  names(records) <- c("year", "prob")

  # Fit the model function using the adapted code from Brook et al. 2019
  model_output <- bbj.2018(
    dd = records,
    iter = n.iter,
    ey = test.time,
    m = "ole",
    plot = FALSE,
    alpha = alpha,
    cores = cores
  )

  # Rename columns back to original
  names(records) <- c("time", "certainty")

  output <- list(
    records = records,
    alpha = alpha,
    test.time = test.time,
    p.extant = as.numeric(model_output$end_year),
    estimate = model_output$MTE,
    conf.int = c(model_output$LCI, model_output$UCI)
  )

  return(output)
}
