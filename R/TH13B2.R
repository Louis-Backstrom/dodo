#' @title Thompson et al.'s (2013) "Quenched" model
#'
#' @description
#' The quenched variant (\eqn{P_Q(X_T|\bm{s})}) of Equation 9 in Thompson et
#' al. 2013. Estimates a posterior probability that the species is extant at
#' the end of the observation period.
#'
#' @param records sighting records in `umcb` format (see
#' \code{\link{convert_dodo}} for details).
#' @param priors nested `list` with two elements: `p` and `q`, each of which is
#' a `list` with the same number of elements as there are sighting classes. For
#' the quenched model, each element is the a two-element vector with the lower
#' and upper bounds of a uniform prior for the rate \eqn{p_i^\alpha} or
#' \eqn{q_i^\alpha}.
#' @param certain which sighting types are considered certain. Defaults to
#' `c(1)` (i.e. only the first sighting class is certain).
#' @param PXT estimate of \eqn{P(X_T)}. Defaults to `NULL`, in which case
#' \eqn{P(X_T) = T_N / T}.
#' @param PE estimate of \eqn{P_E}. Defaults to `NULL`, in which case
#' \eqn{P_E = 1 / T}.
#' @param n.iter number of iterations to calculate averages over. Defaults to
#' 100,000, which is usually sufficient to ensure accurate estimates.
#' @param pb whether to show a progress bar. Defaults to `FALSE`.
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
#' @seealso [TH13B1()]; [TH13B3()]
#'
#' @examples
#' # Run the example analysis from Thompson et al. 2013 (Table 1 etc.)
#' TH13B2(
#'   records = thompson_table1,
#'   priors = list(
#'     p = list(p1 = c(0.2, 0.4), p2 = c(0.4, 0.6), p3 = c(0.3, 0.6)),
#'     q = list(q1 = c(1.0, 1.0), q2 = c(0.2, 0.7), q3 = c(0.1, 0.5))
#'   ),
#'   certain = c(1),
#'   n.iter = 1e3
#' )
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' TH13B2(
#'   records = curlew$umcb,
#'   priors = list(
#'     p = list(
#'       p1 = c(0, 1), p2 = c(0, 1), p3 = c(0, 1),
#'       p4 = c(0, 1), p5 = c(0, 1), p6 = c(0, 1)
#'     ),
#'     q = list(
#'       q1 = c(1.0, 1.0), q2 = c(1.0, 1.0), q3 = c(0.0, 1.0),
#'       q4 = c(0.0, 1.0), q5 = c(0.0, 1.0), q6 = c(0.0, 1.0)
#'     )
#'   ),
#'   certain = c(1, 2),
#'   n.iter = 1e3
#' )
#' }
#'
#' @export

TH13B2 <- function(records, priors, certain = 1, PXT = NULL, PE = NULL,
                   n.iter = 1e5, pb = FALSE) {
  # Sort records
  records <- sort_by(records, ~time)

  # Remove time column and format as matrix
  records_matrix <- as.matrix(subset(records, select = -time))

  # Determine bigT and nA (number of record classes)
  bigT <- nrow(records_matrix)
  nA <- ncol(records_matrix)

  # Determine TN
  if (length(certain) == 1) {
    TN <- max(which(records_matrix[, certain] == 1))
  } else {
    TN <- max(which(rowMaxs(records_matrix[, certain, drop = FALSE]) == 1))
  }

  # Calculate P(X_T) and P(E)
  if (is.null(PXT)) {
    PXT <- TN / bigT
  }

  if (is.null(PE)) {
    PE <- 1 / bigT
  }

  log_PXT <- log(PXT)
  log_PE <- log(PE)

  # Precompute prior bounds
  p_lower <- numeric(nA)
  p_range <- numeric(nA)
  q_lower <- numeric(nA)
  q_range <- numeric(nA)

  for (a in seq_len(nA)) {
    p_lower[a] <- priors$p[[a]][1]
    p_range[a] <- priors$p[[a]][2] - p_lower[a]
    q_lower[a] <- priors$q[[a]][1]
    q_range[a] <- priors$q[[a]][2] - q_lower[a]
  }

  # Run quenched loop
  values <- numeric(n.iter)

  if (pb == TRUE) {
    pbar <- txtProgressBar(min = 0, max = n.iter, style = 3, width = 50)
  }

  for (run in seq_len(n.iter)) {
    log_p_vector <- numeric(bigT)
    log_q_vector <- numeric(bigT)

    # Loop over attributes
    for (a in seq_len(nA)) {
      col_vals <- records_matrix[, a]

      # Create p vector
      u <- runif(bigT)
      vals <- p_lower[a] + p_range[a] * u

      idx <- (col_vals == 1)
      if (any(idx)) {
        vals[idx] <- 1 - vals[idx]
      }

      log_p_vector <- log_p_vector + log(vals)

      # Create q vector
      u <- runif(bigT)
      vals <- q_lower[a] + q_range[a] * u

      if (any(idx)) {
        vals[idx] <- 1 - vals[idx]
      }

      log_q_vector <- log_q_vector + log(vals)
    }

    # Calculate cumulative sums
    cumsum_p <- cumsum(log_p_vector)
    rev_cumsum_q <- rev(cumsum(rev(log_q_vector)))

    # Calculate numerator in log-space
    log_numerator <- cumsum_p[bigT] + log_PXT

    # Calculate denominator in log-space
    if (TN >= bigT) {
      log_denominator <- log_numerator
    } else {
      js <- (TN + 1):bigT

      log_terms <- cumsum_p[js - 1] + rev_cumsum_q[js] + log_PE

      m <- max(log_numerator, log_terms)

      log_denominator <- m + log(exp(log_numerator - m) +
        sum(exp(log_terms - m)))
    }

    # Calculate P_Q(X_T|s) from Equation 9
    values[run] <- exp(log_numerator - log_denominator)

    if (pb == TRUE) {
      setTxtProgressBar(pbar, run)
    }
  }

  if (pb == TRUE) {
    close(pbar)
  }

  # Calculate quenched average
  p.extant <- mean(values)

  # Output
  output <- list(
    records = records,
    priors = priors,
    certain = certain,
    PXT = PXT,
    PE = PE,
    n.iter = n.iter,
    pb = pb,
    p.extant = p.extant
  )

  return(output)
}

# Declare time as a known global variable (column in records data.frame):
utils::globalVariables("time")
