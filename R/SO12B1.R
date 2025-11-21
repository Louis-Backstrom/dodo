#' @title Solow et al.'s (2012) "Uncertain" model
#'
#' @description
#' The model from Solow et al. 2012. Estimates a Bayes factor comparing
#' competing hypotheses of extinction / persistence.
#'
#' @param records sighting records in `ubin` format (see
#' \code{\link{convert_dodo}} for details). Note that all uncertain sightings
#' must follow the final certain sighting.
#' @param init.time start of the observation period.
#' @param increment step size used for integration. Defaults to 0.01, following
#' the original code from Solow & Beet 2014.
#' @param gamma parameter for the exponential prior. Defaults to 1.
#'
#' @returns a `list` object with the original parameters and the Bayes factor
#' included as elements.
#'
#' @note
#' Sampling effort is assumed to be constant. Uses an exponential prior.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A., Smith, W., Burgman, M., Rout, T., Wintle, B., & Roberts, D.
#' (2012). Uncertain sightings and the extinction of the Ivory-billed
#' Woodpecker. *Conservation Biology*, 26(1), 180-184.
#' \doi{10.1111/j.1523-1739.2011.01743.x}
#'
#' **Other References**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#'
#' @examples
#' # Run the Ivory-billed Woodpecker analysis from Solow et al. 2012
#' # Create modified dataset
#' woodpecker$ubin$record <- woodpecker$ubin$certain +
#'   woodpecker$ubin$uncertain
#' woodpecker$ubin$certain <- ifelse(seq_len(nrow(woodpecker$ubin)) <= 48,
#'   woodpecker$ubin$record, 0
#' )
#' woodpecker$ubin$uncertain <- ifelse(seq_len(nrow(woodpecker$ubin)) > 48,
#'   woodpecker$ubin$record, 0
#' )
#' woodpecker$ubin$record <- NULL
#' # Run analysis
#' SO12B1(records = woodpecker$ubin, init.time = 1897, increment = 0.01)
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' SO12B1(curlew$ubin, init.time = 1817, increment = 0.01)
#' }
#'
#' @export

SO12B1 <- function(records, init.time, increment = 0.01, gamma = 1) {
  # Set up model parameters
  certain <- which(records$certain == 1)
  uncertain <- which(records$uncertain == 1)

  DATA <- list()
  DATA$data <- data.frame(
    Year = c(certain, uncertain) + init.time - 1,
    Sightings = c(rep(1, length(certain)), rep(2, length(uncertain)))
  )
  DATA$data <- sort_by(DATA$data, ~Year)

  inputs <- list(T = init.time + nrow(records) - 1, posteriorPrior = "Exp")

  # Fit the model
  fit <- sb14.extended.model(
    DATA = DATA, inputs = inputs, modelnumber = 1,
    increment = increment, increment2 = increment,
    gamma = gamma, SO12 = TRUE
  )

  # Output
  output <- list(
    records = records,
    init.time = init.time,
    increment = increment,
    gamma = gamma,
    Bayes.factor = fit$Results
  )

  return(output)
}
