#' @title Burgman et al.'s (1995) "Discrete-time" model
#'
#' @description
#' Equation 2 from Burgman et al. 1995. Estimates a p-value for
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
#' @seealso [BU95F2()], [BU95F3()]
#'
#' @examples
#' # Run the example analysis from Burgman 1995 (Figure 1b)
#' BU95F1(burgman_figure1b)
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' BU95F1(curlew$cdis)
#' }
#'
#' @export

BU95F1 <- function(records) {
  # Determine number of records
  N <- sum(records)

  # Determine the period in which the last record occurred
  Ce <- max(which(records > 0))

  # Determine the total number of sighting intervals
  CT <- length(records)

  # Calculate p-value
  p.value <- (Ce / CT)^N

  # Output
  output <- list(
    records = records,
    p.value = p.value
  )

  return(output)
}
