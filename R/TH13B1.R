#' @title Thompson et al.'s (2013) "Deterministic" model
#'
#' @description
#' The deterministic variant (\eqn{P_D(X_T|\bm{s})}) of Equation 9 in Thompson
#' et al. 2013. Estimates a posterior probability that the species is extant at
#' the end of the observation period.
#'
#' @param records sighting records in `umcd` format (see
#' \code{\link{convert_dodo}} for details).
#' @param priors nested `list` with two elements: `p` and `q`, each of which is
#' a `list` with the same number of elements as there are sighting classes. For
#' this deterministic model, each element is the point prior for the rate
#' \eqn{p_i^\alpha} or \eqn{q_i^\alpha}.
#' @param certain which sighting types are considered certain. Defaults to
#' `c(1)` (i.e. only the first sighting class is certain).
#' @param PXT estimate of \eqn{P(X_T)}. Defaults to `NULL`, in which case
#' \eqn{P(X_T) = TN / T}.
#' @param PE estimate of \eqn{P_E}. Defaults to `NULL`, in which case
#' \eqn{P_E = 1 / T}.
#'
#' @returns a `list` object with the original parameters and p(extant) included
#' as elements.
#'
#' @note
#' Sampling effort is assumed to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Thompson, C. J., Lee, T. E., Stone, L., McCarthy, M. A., & Burgman, M. A.
#' (2013). Inferring extinction risks from sighting records.
#' *Journal of Theoretical Biology*, 338, 16-22.
#' \doi{10.1016/j.jtbi.2013.08.023}
#'
#' @seealso [TH13B2()]; [TH13B3()]
#'
#' @examples
#' # Run the example analysis from Thompson et al. 2013 (Table 1 etc.)
#' TH13B1(
#'   records = thompson_table1,
#'   priors = list(
#'     p = list(p1 = 0.30, p2 = 0.50, p3 = 0.45),
#'     q = list(q1 = 1.00, q2 = 0.45, q3 = 0.30)
#'   ),
#'   certain = c(1)
#' )
#' # Run an example analysis using the Slender-billed Curlew data
#' \dontrun{
#' TH13B1(
#'   records = curlew$umcd,
#'   priors = list(
#'     p = list(p1 = 0.5, p2 = 0.5, p3 = 0.5, p4 = 0.5, p5 = 0.5, p6 = 0.5),
#'     q = list(q1 = 1.0, q2 = 1.0, q3 = 0.5, q4 = 0.5, q5 = 0.5, q6 = 0.5)
#'   ),
#'   certain = c(1, 2)
#' )
#' }
#'
#' @export

TH13B1 <- function(records, priors, certain = c(1), PXT = NULL, PE = NULL) {
  # Sort records
  records <- sort_by(records, ~time)

  # Remove time column and format as matrix
  records_matrix <- as.matrix(subset(records, select = -time))

  # Determine TN
  if (length(certain) == 1) {
    TN <- max(which(as.numeric(records_matrix[, certain]) == 1))
  } else {
    TN <- max(which(as.numeric(apply(records_matrix[, certain], 1, max)) == 1))
  }

  # Determine bigT
  bigT <- nrow(records_matrix)

  # Calculate P(X_T) and P(E)
  if (is.null(PXT)) {
    PXT <- TN / bigT
  }

  if (is.null(PE)) {
    PE <- 1 / bigT
  }

  # Create p matrix and vector
  p_matrix <- matrix(nrow = nrow(records_matrix), ncol = ncol(records_matrix))

  for (i in 1:nrow(p_matrix)) {
    for (a in 1:ncol(p_matrix)) {
      if (records_matrix[i, a] == 1) {
        p_matrix[i, a] <- 1 - priors$p[[a]]
      } else {
        p_matrix[i, a] <- priors$p[[a]]
      }
    }
  }

  p_vector <- apply(p_matrix, 1, prod)

  # Create q matrix and vector
  q_matrix <- matrix(nrow = nrow(records_matrix), ncol = ncol(records_matrix))

  for (i in 1:nrow(q_matrix)) {
    for (a in 1:ncol(q_matrix)) {
      if (records_matrix[i, a] == 1) {
        q_matrix[i, a] <- 1 - priors$q[[a]]
      } else {
        q_matrix[i, a] <- priors$q[[a]]
      }
    }
  }

  q_vector <- apply(q_matrix, 1, prod)

  # Turn into mpfr objects
  p_vector <- Rmpfr::mpfr(p_vector, precBits = 64)
  q_vector <- Rmpfr::mpfr(q_vector, precBits = 64)

  # Calculate components of Equation 9
  numerator <- prod(p_vector) * PXT

  sums <- list()
  for (j in (TN + 1):bigT) {
    sums[[j]] <- prod(p_vector[1:(j - 1)]) * prod(q_vector[j:bigT]) * PE
  }
  sums <- Filter(Negate(is.null), sums)
  denominator <- numerator + sum(do.call(c, sums))

  # Calculate P_D(X_T|s) from Equation 9
  p.extant <- as.numeric(numerator / denominator)

  # Output
  output <- list(
    records = records,
    priors = priors,
    certain = certain,
    p.extant = p.extant
  )

  return(output)
}

# Declare time as a known global variable (column in records data.frame):
utils::globalVariables("time")
