#' @title Strauss & Sadler's (1989) "Bayesian" model
#'
#' @description
#' Equation 26 from Strauss & Sadler 1989, assuming the prior distribution
#' from equation 22. Estimates a posterior distribution on time of extinction,
#' with associated point estimate and one-sided credible interval.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param t.max maximum time to estimate posterior density out to (defaults to
#' the square of most recent record time; this may not be appropriate for all
#' datasets).
#' @param length.out number of posterior samples to generate (defaults to 10
#' million).
#'
#' @returns a `list` object with the original parameters and the point estimate
#' and credible interval included as elements. The credible interval is a
#' two-element numeric vector called `cred.int`.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Strauss, D., & Sadler, P. M. (1989). Classical Confidence Intervals and
#' Bayesian Probability Estimates for Ends of Local Taxon Ranges.
#' *Mathematical Geology*, 21(4), 411-421. \doi{10.1007/Bf00897326}
#'
#' @seealso [SS89F1()]
#'
#' @examples
#' # Run an example analysis using the Caribbean Monk Seal data
#' SS89B1(monk_seal)
#'
#' @export

SS89B1 <- function(records, alpha = 0.05, t.max = max(records) ^ 2,
                   length.out = 1e7) {

  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # Calculate u_n
  un <- (max(records) - min(records)) ^ (-n + 2) -
    (1 - min(records)) ^ (-n + 2) - (max(records)) ^ (-n + 2) + 1

  # Set up vector of theta_2 values
  theta <- seq(max(records), t.max, length.out = length.out)

  # Calculate posterior density distribution
  pdf <- (n - 2) * ((theta - min(records)) ^ (-n + 1) -
                      theta ^ (-n + 1)) / un

  # Sample from posterior
  posterior <- sample(theta, prob = pdf, replace = T)

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    t.max = t.max,
    length.out = length.out,
    estimate = mean(posterior),
    cred.int = c(max(records), as.numeric(quantile(posterior, 1 - alpha)))
  )

  return(output)

}
