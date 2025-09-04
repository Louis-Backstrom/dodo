#' McInerny et al.'s (2006) "Sighting Rate" model
#'
#' @description
#' Equation 2 from McInerny et al. 2006. Estimates a p-value for testing
#' competing hypotheses of extinction/non-extinction, and a one-tailed
#' \eqn{1 - \alpha} confidence interval on the time of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param remove.first if `init.time` is the time of the first sighting, should
#' that sighting be removed from the record? Defaults to `TRUE`.
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
#' **Other References*
#'
#' Rivadeneira, M. M., Hunt, G., & Roy, K. (2009). The use of sighting records
#' to infer species extinctions: an evaluation of different methods. *Ecology*,
#' 90(5), 1291-1300. \doi{10.1890/08-0316.1}
#'
#' Jaric, I. (2014). The use of sighting records to infer species extinctions:
#' comment. *Ecology*, 95(1), 238. \doi{10.1890/12-1088.1}
#'
#' @export

MC06F1 <- function(records, alpha = 0.05, init.time = min(records),
                   test.time = as.numeric(format(Sys.Date(), "%Y")),
                   remove.first = TRUE) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)
  if (init.time == min(records) & remove.first ==  TRUE) {
    n <- n - 1
  }

  # Calculate p-value
  p.value <- (1 - (n / (max(records) - init.time))) ^
    (test.time - (max(records)))

  # Calculate width of confidence interval
  x <- log(alpha, base = (1 - (n / (max(records) - init.time))))

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.value = p.value,
    conf.int = c(max(records), max(records) + x)
  )

  return(output)

}
