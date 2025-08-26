#' @title Burgman et al.'s (1995) "Runs Test" Model
#'
#' @description
#' Equation 3 from Burgman et al. 1995. Estimates a p-value for
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
#' @export

BU95F2 <- function(records) {

  # Determine the total number of sighting intervals
  CT <- nrow(records)

  # Determine the number of empty cells
  n0 <- sum(records$records == 0)

  # Determine the number of non-empty cells
  n1 <- CT - n0

  # Determine the length of the longest run of empty cells
  r <- max(rle(records$records == 0)$lengths[which(
    rle(records$records == 0)$values == TRUE)])

  # Determine values of k to sum over
  k <- 1:floor(n0/r)

  # Calculate p-value
  part0 <- choose(CT, n0) ^ (-1)
  part1 <- (-1) ^ (k + 1)
  part2 <- choose(n1 + 1, k)
  part3 <- choose(CT - r * k, n1)
  p.value <- part0 * sum(part1 * part2 * part3)
  rm(part0, part1, part2, part3)

  # Output
  output <- list(
    records = records,
    p.value = p.value
  )

  return(output)

}
