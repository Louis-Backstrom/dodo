#' @title Thompson et al.'s (2013) "Annealed" model
#'
#' @description
#' The annealed variant (\eqn{P_A(X_T|\bm{s})}) of Equation 9 in Thompson et
#' al. 2013. Estimates a posterior probability that the species is extant at
#' the end of the survey period.
#'
#' @param records sighting records in `umcd` format (see
#' \code{\link{convert_dodo}} for details).
#' @param priors nested `list` with two elements: `p` and `q`, each of which is
#' a `list` with the same number of elements as there are sighting classes. For
#' the annealed model, each element is the a two-element vector with the lower
#' and upper bounds of a uniform prior for the rate \eqn{p_i^\alpha} or
#' \eqn{q_i^\alpha}.
#' @param certain which sighting types are considered certain. Defaults to
#' `c(1)` (i.e. only the first sighting class is certain).
#' @param PXT estimate of \eqn{P(X_T)}. Defaults to `NULL`, in which case
#' \eqn{P(X_T) = TN / T}.
#' @param PE estimate of \eqn{P_E}. Defaults to `NULL`, in which case
#' \eqn{P_E = 1 / T}.
#' @param n.iter number of iterations to calculate averages over. Defaults to
#' 100,000, which is usually sufficient to ensure accurate estimates.
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
#' @seealso [TH13B1()]; [TH13B2()]
#'
#' @examples
#' # Run the example analysis from Thompson et al. 2013 (Table 1 etc.)
#' TH13B3(
#'   records = thompson_table1,
#'   priors = list(
#'     p = list(p1 = c(0.2, 0.4), p2 = c(0.4, 0.6), p3 = c(0.3, 0.6)),
#'     q = list(q1 = c(1.0, 1.0), q2 = c(0.2, 0.7), q3 = c(0.1, 0.5))
#'   ),
#'   certain = c(1),
#'   n.iter = 1e3
#' )
#' # Run an example analysis using the Slender-billed Curlew data
#' TH13B3(
#'   records = curlew$umcd,
#'   priors = list(
#'     p = list(p1 = c(0, 1), p2 = c(0, 1), p3 = c(0, 1),
#'              p4 = c(0, 1), p5 = c(0, 1), p6 = c(0, 1)),
#'     q = list(q1 = c(1.0, 1.0), q2 = c(1.0, 1.0), q3 = c(0.0, 1.0),
#'              q4 = c(0.0, 1.0), q5 = c(0.0, 1.0), q6 = c(0.0, 1.0))),
#'   certain = c(1, 2),
#'   n.iter = 1e3
#' )
#'
#' @export

TH13B3 <- function(records, priors, certain = 1, PXT = NULL, PE = NULL,
                   n.iter = 1e5) {
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

  # Run annealed loop
  numerator_values <- list()
  denominator_values <- list()

  pb <- txtProgressBar(min = 0, max = n.iter, style = 3, width = 50)
  for (run in 1:n.iter) {
    # Create p matrix and vector
    p_matrix <- matrix(nrow = nrow(records_matrix), ncol = ncol(records_matrix))

    for (i in 1:nrow(p_matrix)) {
      for (a in 1:ncol(p_matrix)) {
        if (records_matrix[i, a] == 1) {
          p_matrix[i, a] <- 1 - runif(n = 1, priors$p[[a]][1], priors$p[[a]][2])
        } else {
          p_matrix[i, a] <- runif(n = 1, priors$p[[a]][1], priors$p[[a]][2])
        }
      }
    }

    p_vector <- apply(p_matrix, 1, prod)

    # Create q matrix and vector
    q_matrix <- matrix(nrow = nrow(records_matrix), ncol = ncol(records_matrix))

    for (i in 1:nrow(q_matrix)) {
      for (a in 1:ncol(q_matrix)) {
        if (records_matrix[i, a] == 1) {
          q_matrix[i, a] <- 1 - runif(n = 1, priors$q[[a]][1], priors$q[[a]][2])
        } else {
          q_matrix[i, a] <- runif(n = 1, priors$q[[a]][1], priors$q[[a]][2])
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

    numerator_values[[run]] <- numerator
    denominator_values[[run]] <- denominator

    setTxtProgressBar(pb, run)
  }
  close(pb)

  # Calculate quenched average
  p.extant <- as.numeric((sum(do.call(c, numerator_values))/n.iter) /
                           (sum(do.call(c, denominator_values))/n.iter))

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
