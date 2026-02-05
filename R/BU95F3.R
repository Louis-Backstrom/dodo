#' @title Burgman et al.'s (1995) "Empty Cells" model
#'
#' @description
#' Equation 4 from Burgman et al. 1995. Estimates a p-value for
#' testing competing hypotheses of extinction/non-extinction.
#'
#' @param records sighting records in `cdis` format (see
#' \code{\link{convert_dodo}} for details).
#'
#' @returns a `list` object with the original parameters and the p-value
#' included as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' Adapts code from the `sExtinct` package
#' (\href{https://github.com/ExperimentalConservation/sExtinct}{GitHub})
#'
#' @references
#' **Key Reference**
#'
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring Threat from
#' Scientific Collections. *Conservation Biology*, 9(4), 923-928.
#' \doi{0.1046/j.1523-1739.1995.09040923.x}
#'
#' **Other References**
#'
#' Clements, C. F., Collen, B., Blackburn, T. M., & Petchey, O. L. (2014).
#' Effects of recent environmental change on accuracy of inferences of
#' extinction status. *Conservation Biology*, 28(4), 971-981.
#' \doi{10.1111/cobi.12329}
#'
#' Grimson, R. C., Aldrich, T. E., & Drane, J. W. (1992). Clustering in sparse
#' data and an analysis of rhabdomyosarcoma incidence.
#' *Statistics in Medicine*, 11(6), 761-768. \doi{10.1002/sim.4780110607}
#'
#' @seealso [BU95F1()], [BU95F2()]
#'
#' @examples
#' # Run the example analysis from Burgman 1995 (Figure 1b)
#' BU95F3(burgman_figure1b)
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' BU95F3(curlew$cdis)
#' }
#'
#' @export

BU95F3 <- function(records) {
  # Heavily adapted from sExtinct (by Christopher Clements) package code!

  # Determine number of records
  N <- sum(records)

  # Determine the total number of sighting intervals
  CT <- length(records)

  # Determine the length of the longest run of empty cells
  r <- max(rle(records == 0)$lengths[which(rle(records == 0)$values == TRUE)])

  # Check if longest run of empty cells occurs at the end, and warn if not
  z <- rle(records == 0)

  z_check <- tail(z$values, 1) && tail(z$lengths, 1) == r

  if (z_check == FALSE) {
    warning("Longest run of empty cells does not occur at the end!")
  }

  mp1 <- Rmpfr::mpfr(1, precBits = 64)
  mp0 <- Rmpfr::mpfr(0, precBits = 64)

  # Calculate sum in for loop (from sExtinct code)
  full_sum <- mp0

  for (j in 1:N) {
    TT <- 0:(j - 1)
    KK <- 0:j

    dummy <- mp0
    for (indexi in 1:length(KK)) {
      i <- indexi - 1
      dummy <- dummy + (((-1)^i) * Rmpfr::chooseMpfr(j, i) *
        ((j - Rmpfr::mpfr(i, precBits = 64))^N))
    }
    Sn <- dummy / factorial(j)

    for (k in 1:(j + 1)) {
      if (k <= (CT / r)) {
        ff <- mp1
        for (tt in TT) {
          term <- CT - (r * k) - tt
          if (term <= 0) {
            ff <- mp0
            break
          }
          ff <- ff * term
        }
        calc <- ((-1)^(k + 1)) * choose(j + 1, k) * ff * Sn
      } else {
        calc <- mp0
      }
      full_sum <- full_sum + calc
    }
  }

  # Calculate p-value
  p.value <- as.numeric(Rmpfr::mpfr(CT, precBits = 64)^-N * full_sum)

  # Output
  output <- list(
    records = records,
    p.value = p.value
  )

  return(output)
}
