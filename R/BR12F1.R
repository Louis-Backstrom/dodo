#' Bradshaw et al.'s (2012) "Inverse-weighted McInerny" model
#'
#' @description
#' Equation 4 from Bradshaw et al. 2012. Estimates a point estimate on the time
#' of extinction.
#'
#' @param records sighting records in `ccon` format (see
#' \code{\link{convert_dodo}} for details).
#' @param alpha desired significance level (defaults to \eqn{\alpha = 0.05}) to
#' use in the underlying McInerny et al. 2006 model.
#'
#' @returns a `list` object with the original parameters and the point estimate
#' included as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant. This is *not* the Gaussian-resampled model (GRIWM) also
#' presented in Bradshaw et al. 2012 (which incorporates radiometric dating
#' error). Duplicate sightings in a single time unit are discarded.
#'
#' @references
#' **Key Reference**
#'
#' Bradshaw, C. J. A., Cooper, A., Turney, C. S. M., & Brook, B. W. (2012).
#' Robust estimates of extinction time in the geological record. *Quaternary*
#' *Science Reviews*, 33, 14-19. \doi{10.1016/j.quascirev.2011.11.021}
#'
#' **Other References**
#'
#' McInerny, G. J., Roberts, D. L., Davy, A. J., & Cribb, P. J. (2006).
#' Significance of sighting rate in inferring extinction and threat.
#' *Conservation Biology*, 20(2), 562-567.
#' \doi{10.1111/j.1523-1739.2006.00377.x}
#'
#' @examples
#' # Run the Woolly Mammoth analysis from Bradshaw et al. 2012
#' BR12F1(mammoth)
#' # Run an example analysis using the Slender-billed Curlew data
#' BR12F1(curlew$ccon)
#'
#' @export

BR12F1 <- function(records, alpha = 0.05) {
  # Sort and de-duplicate records
  records <- sort(unique(records), decreasing = TRUE)

  # Determine number of records
  n <- length(records)

  # Calculate thetas
  theta_k <- c()
  for (i in 2:n) {
    theta_k[i - 1] <- suppressWarnings(
      log(alpha) / log(1 - (i / (max(records[1:i]) - min(records[1:i]))))
    )
  }

  # Replace NA theta values with 0. NA values occur when the terminal dates in
  # the sighting record occur in successive years, in which case the McInerny
  # estimate is (strictly speaking) undefined, but can be safely inferred to be
  # 0 (i.e. extinction immediately following the final record).

  if (sum(is.na(theta_k)) > 0) {
    warning(paste0("Replacing ", sum(is.na(theta_k)),
                   " NA theta_k values with 0."))
  }
  theta_k[is.na(theta_k)] <- 0

  # Calculate weights
  d <- 1 / (records[1] - records[2:n])
  d <- d / max(d)

  # Calculate point estimate
  estimate <- sum(theta_k * d) / sum(d)

  # Output
  output <- list(
    records = sort(records),
    alpha = alpha,
    estimate = estimate + max(records)
  )

  return(output)
}
