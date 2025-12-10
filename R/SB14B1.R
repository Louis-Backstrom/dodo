#' @title Solow & Beet's (2014) "Model 1" model
#'
#' @description
#' The first model from Solow & Beet 2014. Estimates a Bayes factor comparing
#' competing hypotheses of extinction / persistence.
#'
#' @param records sighting records in `ubin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param init.time start of the observation period.
#' @param increment step size used for integration. Defaults to 0.01, following
#' the original code from Solow & Beet 2014.
#'
#' @returns a `list` object with the original parameters and the Bayes factor
#' included as elements.
#'
#' @note
#' Sampling effort is assumed to be constant. Uses a uniform prior.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Beet, A. R. (2014). On uncertain sightings and inference
#' about extinction. *Conservation Biology*, 28(4), 1119-1123.
#' \doi{10.1111/cobi.12309}
#'
#' **Other References**
#'
#' Solow, A., Smith, W., Burgman, M., Rout, T., Wintle, B., & Roberts, D.
#' (2012). Uncertain sightings and the extinction of the Ivory-billed
#' Woodpecker. *Conservation Biology*, 26(1), 180-184.
#' \doi{10.1111/j.1523-1739.2011.01743.x}
#'
#' Carlson, C. J., Bond, A. L., & Burgio, K. R. (2018). Estimating the
#' extinction date of the thylacine with mixed certainty data.
#' *Conservation Biology*, 32(2), 477-483. \doi{10.1111/cobi.13037}
#'
#' Carlson, C. J., Burgio, K. R., Dallas, T. A., & Bond, A. L. (2018). Spatial
#' extinction date estimation: a novel method for reconstructing spatiotemporal
#' patterns of extinction and identifying potential zones of rediscovery.
#' *bioRxiv preprint*. \doi{10.1101/279679}
#'
#' Kodikara, S., Demirhan, H., & Stone, L. (2018). Inferring about the
#' extinction of a species using certain and uncertain sightings.
#' *Journal of Theoretical Biology*, 442, 98-109.
#' \doi{10.1016/j.jtbi.2018.01.015}
#'
#' Burgio, K. R., Carlson, C. J., Bond, A. L., Rubega, M. A., & Tingley, M. W.
#' (2021). The two extinctions of the Carolina Parakeet
#' *Conuropsis carolinensis*. *Bird Conservation International*, 1-8.
#' \doi{10.1017/s0959270921000241}
#'
#' @seealso [SB14B2()]]
#'
#' @examples
#' # Run the Ivory-billed Woodpecker analysis from Solow & Beet 2014
#' SB14B1(records = woodpecker$ubin, init.time = 1897, increment = 0.01)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' SB14B1(curlew$ubin, init.time = 1817, increment = 0.01)
#' }
#'
#' @export

SB14B1 <- function(records, init.time, increment = 0.01) {
  # Set up model parameters
  certain <- which(records$certain == 1)
  uncertain <- which(records$uncertain == 1)

  DATA <- list()
  DATA$data <- data.frame(
    Year = c(certain, uncertain) + init.time - 1,
    Sightings = c(rep(1, length(certain)), rep(2, length(uncertain)))
  )
  DATA$data <- sort_by(DATA$data, ~Year)

  inputs <- list(T = init.time + nrow(records) - 1, posteriorPrior = "Unif")

  # Fit the model
  fit <- sb14.extended.model(
    DATA = DATA, inputs = inputs, modelnumber = 1,
    increment = increment, increment2 = increment,
    gamma = 6, SO12 = FALSE
  ) # NB: gamma parameter not used for Unif!

  # Output
  output <- list(
    records = records,
    init.time = init.time,
    increment = increment,
    Bayes.factor = fit$Results
  )

  return(output)
}
