#' @title McCarthy's (1998) "Partial Solow" Model
#'
#' @description
#' Equation 2 from McCarthy 1998. Estimates a p-value for testing
#' competing hypotheses of extinction/non-extinction.
#'
#' @param records `data.frame` with three columns: `time`, `records`, and
#' `effort`. The `time` column must extend from the start of the observation
#' period (which may be prior to the first sighting) to the end (typically the
#' present day), with evenly-spaced temporal intervals (typically years).
#'
#' @returns a `list` object with the original parameters and the p-value
#' included as elements.
#'
#' @note
#' All sighting records are assumed to be certain.
#'
#' @references
#' **Key Reference**
#'
#' McCarthy, M. A. (1998). Identifying declining and threatened species with
#' museum data. *Biological Conservation*, 83(1), 9-17.
#' \doi{10.1016/S0006-3207(97)00048-7}
#'
#' @export

MC98F1 <- function(records) {

  # Determine number of records
  N <- sum(records$records)

  # Determine the period in which the last record occurred
  t <- max(which(records$records > 0))

  # Determine the pre-terminal collection effort
  et <- sum(records[1:t, ]$effort)

  # Determine the total collection effort
  eT <- sum(records$effort)

  # Calculate p-value
  p.value <- (et / eT) ^ N

  # Output
  output <- list(
    records = records,
    p.value = p.value
  )

  return(output)

}
