#' @title Marshall's (1997) "Non-random" model
#'
#' @description
#' Equations 3 to 5 from Marshall 1997. Estimates a one-sided \eqn{1 - \alpha}
#' confidence interval on time of extinction.
#'
#' @param records sighting records in `cdis` format (see
#' \code{\link{convert_dodo}} for details).
#' @param effort a vector of effort data, of the same length as `records`.
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} confidence interval.
#' @param init.time start of the observation period.
#'
#' @returns a `list` object with the original parameters and the confidence
#' interval included as elements. The confidence interval is a two-element
#' numeric vector called `conf.int`.
#'
#' @note
#' All sighting records are assumed to be certain. Uses the discrete form as
#' presented in Rivadeneira et al. (2009).
#'
#' @references
#' **Key Reference**
#'
#' Marshall, C. R. (1997). Confidence intervals on stratigraphic ranges with
#' nonrandom distributions of fossil horizons. *Paleobiology*, 23(2), 165-173.
#' \doi{10.1017/S0094837300016766}
#'
#' **Other References**
#'
#' Rivadeneira, M. M., Hunt, G., & Roy, K. (2009). The use of sighting records
#' to infer species extinctions: an evaluation of different methods. *Ecology*,
#' 90(5), 1291-1300. \doi{10.1890/08-0316.1}
#'
#' @examples
#' # Run an example analysis using the Lord Howe Gerygone data
#' MA97F1(records = gerygone, effort = gerygone_effort, init.time = 1788)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' MA97F1(records = curlew$cdis, effort = curlew_effort, init.time = 1817)
#' }
#' @export

MA97F1 <- function(records, effort, alpha = 0.05, init.time) {
  # Determine number of records
  N <- sum(records)

  lambda <- alpha^(-1 / (N - 1)) - 1

  # Calculate total effort prior to last sighting
  sum <- sum(effort[1:max(which(records > 0))])
  lsum <- (1 + lambda) * sum

  # Determine the bounding times of the upper confidence interval bound
  lower <- max(which(cumsum(effort) < lsum))
  upper <- min(which(cumsum(effort) > lsum))

  # Determine confidence interval bounds
  conf.int.lower <- init.time + max(which(records > 0)) - 1
  fraction <- (lsum - cumsum(effort)[lower]) /
    (cumsum(effort)[upper] - cumsum(effort)[lower])
  conf.int.upper <- init.time + lower + fraction - 1

  # If upper CI bound is NA (i.e. beyond end of observation period), return Inf
  conf.int.upper <- ifelse(is.na(conf.int.upper), Inf, conf.int.upper)

  # Output
  output <- list(
    records = records,
    effort = effort,
    alpha = alpha,
    init.time = init.time,
    conf.int = c(conf.int.lower, conf.int.upper)
  )

  return(output)
}
