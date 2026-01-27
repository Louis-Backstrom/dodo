#' Thompson et al.'s (2019) "Bayesian Updating" model
#'
#' @description
#' The model from Thompson et al. 2019. Estimates a probability that the
#' species is extant at the test time, and a cumulative Bayes Factor at that
#' point.
#'
#' @param records sighting records in `iucn` format (see
#' \code{\link{convert_dodo}} for details).
#' @param surveys survey effort data (incorporating both dedicated and passive
#' survey effort). A `data.frame` with 11 columns: `time`, `survey`, `epsilon`,
#' `epsilon_lower`, `epsilon_upper`, `p_i`, `p_i_lower`, `p_i_upper`, `p_r`,
#' `p_r_lower`, and `p_r_upper`. The number of rows must match `records`.
#' @param init.time the time point at which to start the Bayesian updating
#' algorithm. At this point, \eqn{P(X_t)} is set to `pi`.
#' @param test.time time point to retrospectively calculate extinction
#' probability at.
#' @param pi prior probability that the species is extant at `init.time`
#' (defaults to \eqn{\pi = 0.5}).
#'
#' @returns a `list` object with the original parameters and the final
#' estimates of p(extant) (\eqn{P(X_t)}) and the cumulative Bayes Factor
#' included as elements.
#'
#' @note
#' This model incorporates both sighting uncertainty and variable survey effort.
#'
#' @references
#' **Key Reference**
#'
#' Thompson, C. J., Kodikara, S., Burgman, M. A., Demirhan, H., & Stone, L.
#' (2019). Bayesian updating to estimate extinction from sequential observation
#' data. *Biological Conservation*, 229, 26-29.
#' \doi{10.1016/j.biocon.2018.11.003}
#'
#' @examples
#' # Run the Alaotra Grebe analysis from Thompson et al. 2019 (Fig. 1)
#' TH19B1(
#'   records = grebe, surveys = grebe_surveys, init.time = 1990,
#'   test.time = 1998, pi = 0.5
#' )
#'
#' @export

TH19B1 <- function(records, surveys, init.time, test.time, pi = 0.5) {
  # Set up merged records and surveys data.frame
  data <- merge(records, surveys, by = "time")
  data$BF <- NA
  data$CBF <- NA
  data$PXt <- NA

  # If test.time is before init.time, return p.extant = 1
  if (test.time < init.time) {
    output <- list(
      records = records,
      surveys = surveys,
      init.time = init.time,
      test.time = test.time,
      pi = pi,
      Bayes.factor = 0,
      p.extant = 1
    )

    return(output)
  }

  # Filter to init.time to test.time
  data <- data[data$time >= init.time, ]
  data <- data[data$time <= test.time, ]

  # Calculate yearly Bayes Factors
  data$BF <- ifelse(data$record == TRUE, 1 - data$certainty,
    (1 - data$epsilon * data$p_i * data$p_r)^-1
  )

  # Calculate cumulative Bayes Factors
  data$CBF <- cumprod(data$BF)

  # Calculate yearly p(extant)
  for (i in 1:nrow(data)) {
    if (i == 1) {
      data$PXt[i] <- pi
    } else {
      data$PXt[i] <- (1 + data$CBF[i - 1] * (1 / pi - 1))^-1
    }
  }

  # Output
  output <- list(
    records = records,
    surveys = surveys,
    init.time = init.time,
    test.time = test.time,
    pi = pi,
    Bayes.factor = tail(data$CBF, 1),
    p.extant = tail(data$PXt, 1)
  )

  return(output)
}
