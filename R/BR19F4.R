#' @title Brook et al.'s (2019) "McInerny" model
#'
#' @description
#' The McInerny et al. 2006 variant of the model from Brook et al. 2019.
#' Estimates a "probability" of extinction at the test time, and a two-tailed
#' \eqn{1 - \alpha} confidence interval and point estimate on the time of
#' extinction. Sighting uncertainty is incorporated.
#'
#' @param records `data.frame` with two columns: `time` and `certainty`. The
#' `time` column contains the date of all sightings, which each have an
#' associated `certainty`, expressed as a probability in the interval \[0, 1]
#' (values close to 1 imply a high likelihood that the observation is correct).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param cores number of cores to use (defaults to `NULL`, in which case
#' `parallel::detectCores()` is run).
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
#' BR19F4(woodpecker3, alpha = 0.1, cores = 2)
#' # NB: alpha = 0.1 as this package presents two-sided confidence intervals
#' # vs. one-sided confidence intervals in the original paper; the upper bound
#' # of the one-sided alpha = 0.05 CI in Brook et al. is equivalent to the
#' # upper bound of the two-sided alpha = 0.1 CI in this example.
#'
#' @export

BR19F4 <- function(records, alpha = 0.05,
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   cores = NULL) {
  # Sort records
  records <- sort_by(records, ~time)

  # Rename columns to fit model function
  names(records) <- c("year", "prob")

  # Fit the model function using the adapted code from Brook et al. 2019
  model_output <- bbj.2018(
    dd = records,
    iter = 1e4,
    ey = test.time,
    m = "mc06",
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
