#' @title Alroy's (2015) "Agnostic" model
#'
#' @description
#' The model from Alroy 2015. Estimates a posterior probability that the
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
#' to be constant. Uses a uniform prior, as suggested by Alroy (pers. comm.).
#'
#' @references
#' **Key Reference**
#'
#' Alroy, J. (2015). Current extinction rates of reptiles and amphibians.
#' *Proceedings of the National Academy of Sciences of the United States of*
#' *America*, 112(42), 13003-13008. \doi{10.1073/pnas.1508681112}
#'
#' **Other References**
#'
#' Alroy, J. (2016). On a conservative Bayesian method of inferring extinction.
#' *Paleobiology*, 42(4), 670-679. \doi{10.1017/pab.2016.12}
#'
#' Alroy, J. (2016). A simple Bayesian method of inferring extinction: reply.
#' *Ecology*, 97(3), 798-800. \doi{10.1002/ecy.1321}
#'
#' @seealso [AL14B1()]]
#'
#' @examples
#' # Run an example analysis using the Caribbean Monk Seal data
#' AL15B1(as.integer(1915:1992 %in% monk_seal))
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' AL15B1(curlew$cbin)
#' }
#'
#' @export

AL15B1 <- function(records) {
  # Run agnostic model function
  model_output <- agnostic(records, uniform = TRUE) # uniform = TRUE per Alroy

  # Output
  output <- list(
    records = records,
    p.extant = 1 - tail(model_output, 1)
  )

  return(output)
}

#' @title Agnostic model function from Alroy (2015).
#'
#' @description
#' Helper function. From provided code (agnostic.R)
#'
#' @param onoff the 1/0 sighting vector (i.e. records in `cbin` format).
#' @param total not used.
#' @param prior not used.
#' @param uniform whether or not to use a uniform prior.
#'
#' @returns a `numeric` vector of poster probabilities of extinction for each
#' discrete time interval in the survey per.
#'
#' @references
#' **Key Reference**
#'
#' Alroy, J. (2015). Current extinction rates of reptiles and amphibians.
#' *Proceedings of the National Academy of Sciences of the United States of*
#' *America*, 112(42), 13003-13008. \doi{10.1073/pnas.1508681112}
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

agnostic <- function(onoff, total, prior, uniform) {
  # Intelligently impute missing parameter values
  if (missing(total) || total == "") {
    total <- rep(1, length(onoff))
  }
  if (sum(total) == 0) {
    return(rep(NA, length(onoff)))
  }
  if (missing(prior) || prior == "") {
    prior <- 0.5
  }
  if (missing(uniform) || uniform != T) {
    uniform <- FALSE
  }

  totalBins <- length(onoff)
  bins <- length(onoff)
  for (i in 1:bins) {
    if (total[i] > 0) {
      lastControl <- i
    }
  }
  bins <- lastControl
  onoff <- onoff[1:bins]
  total <- total[1:bins]
  if (onoff[bins] > 0) {
    return(rep(0, totalBins))
  }

  withData <- 0
  for (i in 1:bins) {
    if (onoff[i] > 0) {
      withData <- withData + 1
    }
  }
  if (withData <= 1) {
    return(rep(NA, totalBins))
  }

  first <- 0
  last <- 0
  for (i in 1:bins) {
    if (onoff[i] > 0) {
      if (first == 0) {
        first <- i
      }
      last <- i
    }
  }

  n <- sum(onoff)

  # Set up and calculate likelihoods
  c <- choose(sum(total[1:last]) - 1, n - 1)
  aliveLike <- c / choose(sum(total), n)

  extinctLike <- rep(0, totalBins)
  mu <- log(1 - prior) / bins
  for (i in (last + 1):bins) {
    if (uniform == TRUE) {
      extinctLike[i] <- prior * c / choose(sum(total[1:(i - 1)]), n) / bins
    } else {
      extinctLike[i] <- c / choose(sum(total[1:(i - 1)]), n) *
        (exp((i - 1) * mu) - exp(i * mu))
    }
  }

  # Calculate posteriors from likelihoods and prior
  post <- extinctLike / (sum(extinctLike) + (1 - prior) * aliveLike)

  post <- cumsum(post)
  if (bins < totalBins) {
    post[(bins + 1):totalBins] <- post[length(post)]
  }

  return(post)
}
