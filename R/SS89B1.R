#' @title Strauss & Sadler's (1989) "Bayesian" model
#'
#' @description
#' Equation 26 from Strauss & Sadler 1989, assuming the prior distribution
#' from equation 22. Estimates a posterior distribution on time of extinction,
#' with associated point estimate and one-sided credible interval.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param length.out number of posterior samples to generate (defaults to 100
#' thousand).
#' @param scale factor to scale sighting records by. Defaults to 0.01; adjust
#' if warned.
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
#' SS89B1(monk_seal, length.out = 1e5)
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' SS89B1(curlew$ccon, length.out = 1e5)
#' }
#'
#' @export

SS89B1 <- function(records, alpha = 0.05, length.out = 1e7, scale = 0.01) {
  # Sort records
  records <- sort(records)

  # Scale records
  records_scaled <- scale + scale * (records - min(records)) /
    (max(records) - min(records))

  # Determine number of records
  n <- length(records_scaled)

  # Calculate u_n
  parta <- Rmpfr::mpfr(max(records_scaled) - min(records_scaled),
    precBits = 1024
  )^(-n + 2)
  partb <- (1 - Rmpfr::mpfr(min(records_scaled), precBits = 1024))^(-n + 2)
  partc <- (Rmpfr::mpfr(max(records_scaled), precBits = 1024))^(-n + 2)
  un <- parta - partb - partc + 1

  # Set up vector of theta_2 values
  theta <- Rmpfr::mpfr(seq(max(records_scaled), 1, length.out = length.out),
    precBits = 1024
  )

  # Calculate posterior density distribution
  pdf <- (n - 2) * ((theta - min(records_scaled))^(-n + 1) -
    theta^(-n + 1)) / un

  # Sample from posterior
  posterior <- sample(as.numeric(theta), prob = as.numeric(pdf), replace = T)

  # Check if close to scale limit
  if (max(posterior) > 0.9) {
    warning(paste0(
      "Scale factor may be too large - try setting lower (e.g. ",
      scale * 0.1, ")"
    ))
  }

  # Back-transform estimates
  estimate <- (mean(posterior) - scale) * ((max(records) - min(records)) /
    scale) + min(records)
  cred.int.upper <- (as.numeric(quantile(posterior, 1 - alpha)) - scale) *
    ((max(records) - min(records)) / scale) + min(records)

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    length.out = length.out,
    scale = scale,
    estimate = estimate,
    cred.int = c(max(records), cred.int.upper)
  )

  gc()
  return(output)
}
