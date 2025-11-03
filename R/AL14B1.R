#' @title Alroy's (2014) "Creeping Shadow of a Doubt" model
#'
#' @description
#' The model from Alroy 2014. Estimates a posterior probability that the
#' species is extant at the end of the observation period.
#'
#' @param records sighting records in `cbin` format (see
#' \code{\link{convert_dodo}} for details).
#'
#' @returns a `list` object with the original parameters and p(extant) included
#' as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Alroy, J. (2014). A simple Bayesian method of inferring extinction.
#' *Paleobiology*, 40(4), 584-607. \doi{10.1666/13074}
#'
#' **Other References**
#'
#' Alroy, J. (2016). On a conservative Bayesian method of inferring extinction.
#' *Paleobiology*, 42(4), 670-679. \doi{10.1017/pab.2016.12}
#'
#' Alroy, J. (2016). A simple Bayesian method of inferring extinction: reply.
#' *Ecology*, 97(3), 798-800. \doi{10.1002/ecy.1321}
#'
#' @seealso [AL15B1()]]
#'
#' @examples
#' # Run an example analysis using the Caribbean Monk Seal data
#' AL14B1(as.integer(1915:1992 %in% monk_seal))
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' AL14B1(curlew$cbin)
#' }
#'
#' @export

AL14B1 <- function(records) {
  # Convert records to matrix format
  pa <- as.matrix(t(records))

  # Specify model parameters
  years <- ncol(pa)
  n <- 1

  abs <- 0
  pres <- 0
  lastseen <- 0
  post <- rep(0, years)

  # Run model
  for (y in 1:years) {
    if (pa[n, y] > 0) {
      embedded <- abs
      pres <- pres + 1
      post[y] <- 0
      lastseen <- y
    } else {
      abs <- abs + 1
    }

    if (pres > 2 && lastseen < y) {
      qA <- max(embedded, 1) / (y - 1)
      # Strauss and Sadler + exponential prior
      r <- 2 * lastseen * (pres + 1) / (pres - 1)
      pE0 <- mu(r)
      pE <- 0
      for (i in (lastseen + 1):y) {
        pE <- (pE + (1 - pE) * pE0) /
          (pE + (1 - pE) * pE0 + (1 - pE) * (1 - pE0) * qA)
        post[y] <- pE
      }
    }
  }

  # Output
  output <- list(
    records = records,
    p.extant = 1 - tail(post, 1)
  )

  return(output)
}

#' @title mu function from Alroy (2014).
#'
#' @description
#' Helper function. From provided code (BayesianExtinction_calculation.R)
#'
#' @param range a number.
#'
#' @returns a number.
#'
#' @references
#' **Key Reference**
#'
#' Alroy, J. (2014). A simple Bayesian method of inferring extinction.
#' *Paleobiology*, 40(4), 584-607. \doi{10.1666/13074}
#'
#' **Other References**
#'
#' Alroy, J. (2016). On a conservative Bayesian method of inferring extinction.
#' *Paleobiology*, 42(4), 670-679. \doi{10.1017/pab.2016.12}
#'
#' Alroy, J. (2016). A simple Bayesian method of inferring extinction: reply.
#' *Ecology*, 97(3), 798-800. \doi{10.1002/ecy.1321}
#'
#' @noRd

mu <- function(range) {
  return(1 - exp(log(0.5) / range))
}
