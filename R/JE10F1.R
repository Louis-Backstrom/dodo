#' Jarić & Ebenhard's (2010) "Stationary" model
#'
#' @description
#' Equation 4 from Jarić & Ebenhard 2010. Estimates a p-value for testing
#' competing hypotheses of extinction/non-extinction and a one-sided
#' \eqn{1 - \alpha} confidence interval on the time of extinction.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#'
#' @returns a `list` object with the original parameters and the p-value and
#' confidence interval included as elements. The confidence interval is a
#' two-element numeric vector called `conf.int`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Jarić, I., & Ebenhard, T. (2010). A method for inferring extinction based on
#' sighting records that change in frequency over time. *Wildlife Biology*,
#' 16(3), 267-275. \doi{10.2981/09-044}
#'
#' **Other References**
#'
#' McCrea, R. S., Cheale, T., Campillo-Funollet, E., & Roberts, D. L. (2024).
#' Inferring species extinction from sighting data. *Cambridge Prisms: *
#' *Extinction*, 2, e19. \doi{10.1017/ext.2024.18}
#'
#' @seealso [JE10F2()]
#'
#' @examples
#' # Run the Black-footed Ferret analysis from Jarić & Ebenhard 2010
#' JE10F1(ferret$ccon - 5, test.time = 223) # shift dates to align with paper
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' JE10F1(curlew$ccon, init.time = 1817, test.time = 2022)
#' }
#'
#' @export

JE10F1 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {
  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # If using first record as init.time, remove this from the record sequence
  if (init.time == min(records)) {
    records <- tail(records, -1)
    n <- n - 1
  }

  # Determine length of sighting record
  tn <- max(records) - init.time

  # Determine current time
  bigT <- test.time - init.time

  # Calculate p-value
  p.value <- ((tn - 1) / (n - 1)) /
    (((tn - 1) / (n - 1)) + (bigT - tn))

  # Calculate width of confidence interval
  x <- tn + (((tn - 1) / (n - 1)) + 0) * ((1 - alpha) / alpha)

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.value = p.value,
    conf.int = c(init.time + tn, init.time + x)
  )

  return(output)
}
