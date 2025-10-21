#' McInerny et al.'s (2006) "Sighting Rate" model
#'
#' @description
#' Equation 2 from McInerny et al. 2006. Estimates a p-value for testing
#' competing hypotheses of extinction/non-extinction, and a one-tailed
#' \eqn{1 - \alpha} confidence interval on the time of extinction.
#'
#' @param records sighting records in `cbin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start time for the discrete-time sighting record.
#'
#' @returns a `list` object with the original parameters and the p-value and
#' confidence interval included as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' McInerny, G. J., Roberts, D. L., Davy, A. J., & Cribb, P. J. (2006).
#' Significance of sighting rate in inferring extinction and threat.
#' *Conservation Biology*, 20(2), 562-567.
#' \doi{10.1111/j.1523-1739.2006.00377.x}
#'
#' **Other References**
#'
#' Rivadeneira, M. M., Hunt, G., & Roy, K. (2009). The use of sighting records
#' to infer species extinctions: an evaluation of different methods. *Ecology*,
#' 90(5), 1291-1300. \doi{10.1890/08-0316.1}
#'
#' Jaric, I. (2014). The use of sighting records to infer species extinctions:
#' comment. *Ecology*, 95(1), 238. \doi{10.1890/12-1088.1}
#'
#' @examples
#' # Run an example analysis using the Black-footed Ferret data
#' MC06F1(as.integer(ferret$cdis != 0), init.time = 1972)
#' # Run an example analysis using the Slender-billed Curlew data
#' MC06F1(curlew$cbin, init.time = 1817)
#'
#'
#' @export

MC06F1 <- function(records, alpha = 0.05, init.time) {
  # Determine number of records
  n <- sum(records)

  # Determine time of last record
  tn <- max(which(records == 1))

  # Determine length of sighting record
  bigT <- length(records)

  # Calculate p-value
  p.value <- (1 - (n / tn))^(bigT - tn)

  # Calculate width of confidence interval
  x <- log(alpha, base = (1 - (n / tn)))

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    # test.time = test.time,
    p.value = p.value,
    conf.int = c(init.time + tn - 1, init.time + tn - 1 + x)
  )

  return(output)
}
