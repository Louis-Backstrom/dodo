#' @title Jarić & Roberts' (2014) "Strauss & Sadler 1989" Model
#'
#' @description
#' Modification of the classical confidence interval presented by Strauss &
#' Sadler 1989 from Jarić & Roberts 2014. Estimates a one-tailed
#' \eqn{1 - \alpha} confidence interval on the time of
#' extinction. Sighting uncertainty is incorporated.
#'
#' @param records `data.frame` with two columns: `time` and `certainty`. The
#' `time` column contains the date of all sightings, which each have an
#' associated `certainty`, expressed as a probability in the interval \[0, 1]
#' (values close to 1 imply a high likelihood that the observation is correct).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting.
#'
#' @returns a `list` object with the original parameters and the confidence
#' interval included as elements. The confidence interval is a two-element
#' numeric vector called `conf.int`.
#'
#' @note
#' Sampling effort is assumed to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Jarić, I., & Roberts, D. L. (2014). Accounting for observation reliability
#' when inferring extinction based on sighting records. *Biodiversity and*
#' *Conservation*, 23(11), 2801-2815. \doi{10.1007/s10531-014-0749-8}
#'
#' **Other References**
#'
#' Strauss, D., & Sadler, P. M. (1989). Classical Confidence Intervals and
#' Bayesian Probability Estimates for Ends of Local Taxon Ranges. Mathematical
#' Geology, 21(4), 411-421. \doi{10.1007/Bf00897326}
#'
#' @seealso [JR14F2()], [JR14F3()], [JR14F4()]
#'
#' @export

JR14F1 <- function(records, alpha = 0.05, init.time = min(records$time)) {

  # Sort records
  records <- sort_by(records, ~time)

  # Add relative time column
  records$rtime <- records$time - init.time

  # Calculate r
  r <- sum(records$certainty)

  # Calculate tr
  sum_vector <- c()
  for (i in 1:nrow(records)) {
    sum_vector[i] <- records$rtime[i] *
      records$certainty[i] *
      ifelse(i == nrow(records), 1, prod(1 - records[(i + 1):nrow(records), ]$certainty))
  }
  tr <- sum(sum_vector)

  # Calculate lambda
  lambda <- alpha ^ (-1 / (r - 1)) - 1

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    conf.int = c(init.time + tr, init.time + tr + lambda * tr)
  )

  return(output)

}
