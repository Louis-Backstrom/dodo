#' @title Burgman et al.'s (1995) "Runs Test" model
#'
#' @description
#' Equation 3 from Burgman et al. 1995. Estimates a p-value for
#' testing competing hypotheses of extinction/non-extinction.
#'
#' @param records sighting records in `cdis` format (see
#' \code{\link{convert_dodo}} for details).
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
#' @seealso [BU95F1()], [BU95F3()]
#'
#' @examples
#' # Run the example analysis from Burgman 1995 (Figure 1b)
#' BU95F2(burgman_figure1b)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' BU95F2(curlew$cdis)
#' }
#'
#' @export

BU95F2 <- function(records) {
  # Determine the total number of sighting intervals
  CT <- length(records)

  # Determine the number of empty cells
  n0 <- sum(records == 0)

  # Determine the number of non-empty cells
  n1 <- CT - n0

  # Determine the length of the longest run of empty cells
  r <- max(rle(records == 0)$lengths[which(
    rle(records == 0)$values == TRUE
  )])

  # Check if longest run of empty cells occurs at the end, and warn if not
  z <- rle(records == 0)

  z_check <- tail(z$values, 1) && tail(z$lengths, 1) == r

  if (z_check == FALSE) {
    warning("Longest run of empty cells does not occur at the end!")
  }

  # Determine values of k to sum over
  k <- 1:floor(n0 / r)

  # Calculate p-value
  part0 <- choose(CT, n0)^(-1)
  part1 <- (-1)^(k + 1)
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
