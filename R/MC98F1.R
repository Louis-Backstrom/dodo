#' @title McCarthy's (1998) "Partial Solow" model
#'
#' @description
#' Equation 2 from McCarthy 1998. Estimates a p-value for testing
#' competing hypotheses of extinction/non-extinction.
#'
#' @param records sighting records in `cdis` format (see
#' \code{\link{convert_dodo}} for details).
#' @param effort a vector of effort data, of the same length as `records`.
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
#' @examples
#' # Run an example analysis using the Lord Howe Gerygone data
#' MC98F1(records = gerygone, effort = gerygone_effort)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' MC98F1(records = curlew$cdis, effort = curlew_effort)
#' }
#'
#' @export

MC98F1 <- function(records, effort) {
  # Determine number of records
  N <- sum(records)

  # Determine the period in which the last record occurred
  t <- max(which(records > 0))

  # Determine the pre-terminal collection effort
  et <- sum(effort[1:t])

  # Determine the total collection effort
  eT <- sum(effort)

  # Calculate p-value
  p.value <- (et / eT)^N

  # Output
  output <- list(
    records = records,
    effort = effort,
    p.value = p.value
  )

  return(output)
}
