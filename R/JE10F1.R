#' Jarić & Ebenhard's (2010) "Stationary" model
#'
#' @description
#' Equation 4 from Jarić & Ebenhard 2010. Estimates a p-value for testing
#' competing hypotheses of extinction/non-extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
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
#' Jarić, I., & Ebenhard, T. (2010). A method for inferring extinction based on
#' sighting records that change in frequency over time. *Wildlife Biology*,
#' 16(3), 267-275. \doi{10.2981/09-044}
#'
#' @seealso [JE10F2()]
#'
#' @examples
#' # Run the Black-footed Ferret analysis from Jarić & Ebenhard 2010
#' JE10F1(ferret1 - 5, test.time = 223) # shift dates to align with paper
#'
#' @export

JE10F1 <- function(records, test.time = as.numeric(format(Sys.Date(), "%Y"))) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # Calculate p-value
  p.value <- ((max(records) - 1) / (n - 1)) /
    (((max(records) - 1) / (n - 1)) + (test.time - max(records)))

  # Output
  output <- list(
    records = records,
    test.time = test.time,
    p.value = p.value
  )

  return(output)

}
