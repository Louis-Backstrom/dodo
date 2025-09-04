#' Bradshaw et al.'s (2012) "Inverse-weighted McInerny" model
#'
#' @description
#' Equation 4 from Bradshaw et al. 2012. Estimates a point estimate on the time
#' of extinction.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
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
#' error).
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
#' @export

BR12F1 <- function(records, alpha = 0.05) {

  # Sort records
  records <- sort(records, decreasing = T)

  # Determine number of records
  n <- length(records)

  # Calculate thetas
  theta_k <- c()
  for (i in 2:n) {
    theta_k[i - 1] <- MC06F1(records[1:i], alpha = alpha,
                             remove.first = FALSE)$conf.int[2]
  }

  # Calculate weights
  d <- 1 / (records[1] - records[2:n])
  d <- d / max(d)

  # Calculate point estimate
  estimate <- sum(theta_k * d) / sum(d)

  # Output
  output <- list(
    records = sort(records),
    alpha = alpha,
    estimate = estimate
  )

  return(output)

}
