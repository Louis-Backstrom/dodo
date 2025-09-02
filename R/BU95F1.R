#' @title Burgman et al.'s (1995) "Discrete-time" Model
#'
#' @description
#' Equation 2 from Burgman et al. 1995. Estimates a p-value for
#' testing competing hypotheses of extinction/non-extinction.
#'
#' @param records `data.frame` with two columns: `time` and `records`. The
#' `time` column must extend from the start of the observation period (which
#' may be prior to the first sighting) to the end (typically the present day),
#' with evenly-spaced temporal intervals (typically years).
#'
#' @returns a `list` object with the original parameters and the p-value
#' included as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring Threat from
#' Scientific Collections. *Conservation Biology*, 9(4), 923-928.
#' \doi{0.1046/j.1523-1739.1995.09040923.x}
#'
#' @seealso [BU95F2()], [BU95F3()]
#'
#' @export

BU95F1 <- function(records) {

  # Determine number of records
  N <- sum(records$records)

  # Determine the period in which the last record occurred
  Ce <- max(which(records$records > 0))

  # Determine the total number of sighting intervals
  CT <- nrow(records)

  # Calculate p-value
  p.value <- (Ce / CT) ^ N

  # Output
  output <- list(
    records = records,
    p.value = p.value
  )

  return(output)

}
