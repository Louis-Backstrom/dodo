#' @title Jarić & Roberts' (2014) "Solow 1993a" model
#'
#' @description
#' Equations 7-10 from Jarić & Roberts 2014. Estimates a p-value for testing
#' competing hypotheses of extinction/non-extinction, and a one-tailed
#' \eqn{1 - \alpha} confidence interval and point estimate on the time of
#' extinction. Sighting uncertainty is incorporated.
#'
#' @param records `data.frame` with two columns: `time` and `certainty`. The
#' `time` column contains the date of all sightings, which each have an
#' associated `certainty`, expressed as a probability in the interval \[0, 1]
#' (values close to 1 imply a high likelihood that the observation is correct).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#'
#' @returns a `list` object with the original parameters and the p-value, point
#' estimate, and confidence interval included as elements. The confidence
#' interval is a two-element numeric vector called `conf.int`.
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
#' Solow, A. R. (1993). Inferring Extinction from Sighting Data. *Ecology*,
#' 74(3), 962-964. \doi{10.2307/1940821}
#'
#' @seealso [JR14F1()], [JR14F3()], [JR14F4()]
#'
#' @export

JR14F2 <- function(records, alpha = 0.05, init.time = min(records$time),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {

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

  # Calculate p-value
  p.value <- (tr / (test.time - init.time)) ^ r

  # Calculate point estimate
  estimate <- (tr * (r + 1) / r) + init.time

  # Calculate upper bound of confidence interval
  conf.int.upper <- tr / (alpha ^ (1 / r)) + init.time

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.value = p.value,
    estimate = estimate,
    conf.int = c(init.time + tr, conf.int.upper)
  )

  return(output)

}
