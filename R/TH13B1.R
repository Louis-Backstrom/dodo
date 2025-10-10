#' @title Thompson et al.'s (2013) "Deterministic" model
#'
#' @description
#' The deterministic variant (\eqn{P_D(X_T|\bm{s})}) of Equation 9 in Thompson
#' et al. 2013. Estimates a posterior probability that the species is extant at
#' the end of the survey period.
#'
#' @param records `data.frame` with two or more columns: `time` and one or more
#' columns with the sighting history (0/1) for each distinct sighting class
#' (the column names for these columns is not important).
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

  # Calculate components of Equation 9
  numerator <- prod(p_vector) * PXT

  sums <- c()
  for (j in (TN + 1):bigT) {
    sums[j] <- prod(p_vector[1:(j - 1)]) * prod(q_vector[j:bigT]) * PE
  }
  denominator <- numerator + sum(sums, na.rm = T)

  # Calculate P_D(X_T|s) from Equation 9
  p.extant <- numerator / denominator

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
